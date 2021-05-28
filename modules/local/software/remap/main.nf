
process sam_to_fastq_paired {
	tag "${meta.sampleName}"

	input:
		tuple val (meta), val (type), path (bam)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.R1.fastq.gz"), path ("${meta.sampleName}.${type}.R2.fastq.gz"), emit: result
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.R1.fastq.gz.md5"), path ("${meta.sampleName}.${type}.R2.fastq.gz.md5"), emit: checksum

	script:
	"""#!/usr/bin/env bash

java -Xmx${params.picard.ram}G -Dpicard.useLegacyParser=false -jar ${params.picard.jar} SamToFastq \\
	-INPUT ${bam} \\
	-FASTQ ${meta.sampleName}.${type}.R1.fastq.gz \\
	-SECOND_END_FASTQ ${meta.sampleName}.${type}.R2.fastq.gz \\
	-INCLUDE_NON_PF_READS true \\
	-VALIDATION_STRINGENCY LENIENT

md5sum ${meta.sampleName}.${type}.R1.fastq.gz > ${meta.sampleName}.${type}.R1.fastq.gz.md5
md5sum ${meta.sampleName}.${type}.R2.fastq.gz > ${meta.sampleName}.${type}.R2.fastq.gz.md5

	"""
}

process fastqc_paired {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${meta.sampleName}/results/FastQC", mode: "copy", pattern: "QC/*", saveAs: { it.replaceFirst ("^QC/","") }

	input:
		tuple val (meta), val (type), path (fastq_r1), path (fastq_r2)

	output:
		tuple val (meta), val (type), path (fastq_r1), path (fastq_r2), env (PHRED), emit: result
		path ("QC/*"), emit: fastqc

	script:
	"""#!/usr/bin/env bash
mkdir QC
fastqc -t ${params.fastqc.threads} \\
	${fastq_r1} \\
	${fastq_r2} \\
	--outdir=QC

illumina_version=\$(unzip -p QC/${meta.sampleName}.${type}.R1_fastqc.zip ${meta.sampleName}.${type}.R1_fastqc/fastqc_data.txt | grep Encoding | awk '{print \$5}')

case \${illumina_version} in
	1\\.[8-9])
	PHRED="phred33"
	;;
	*)
	PHRED="phred64"
	;;
esac
echo "PHRED: '\${PHRED}'"
	"""
}

process trim_paired {
	tag "${meta.sampleName}"

	input:
		tuple val (meta), val(type), path (fastq_r1), path (fastq_r2), val (phred)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.R1.passed.fastq.gz"), path ("${meta.sampleName}.${type}.R2.passed.fastq.gz"), emit: result

	script:
	"""#!/usr/bin/env bash

java -Xmx${params.trimmomatic.ram}G -jar ${params.trimmomatic.jar} PE \\
	-threads ${params.trimmomatic.threads} \\
	-${phred} \\
	${fastq_r1} \\
	${fastq_r2} \\
	${meta.sampleName}.${type}.R1.passed.fastq.gz \\
	${meta.sampleName}.${type}.R1.not_passed.fastq.gz \\
	${meta.sampleName}.${type}.R2.passed.fastq.gz \\
	${meta.sampleName}.${type}.R2.not_passed.fastq.gz \\
	LEADING:25 \\
	TRAILING:25 \\
	MINLEN:50 \\
	SLIDINGWINDOW:10:25 \\
	ILLUMINACLIP:${params.trimmomatic.dir}/adapters/TruSeq3-PE-2.fa:2:30:10

	"""
}

process bwa_mem_paired {
	tag "${meta.sampleName}"

	input:
		val (reference_dir)
		tuple val (meta), val (type), path (fastq_r1), path (fastq_r2)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.cleaned.bam"), emit: result

	script:
	"""#!/usr/bin/env bash

bwa mem -t ${params.bwa_mem.threads} \\
	${reference_dir} \\
	-Y \\
	-K ${params.bwa_mem.input_bases} \\
	-v 1 \\
	${fastq_r1} \\
	${fastq_r2} \\
	| java -Xmx${params.picard.ram}G -Dpicard.useLegacyParser=false \\
	-jar ${params.picard.jar} CleanSam \\
	-I /dev/stdin \\
	-O ${meta.sampleName}.${type}.cleaned.bam \\
	-VALIDATION_STRINGENCY LENIENT

	"""
}

process mark_duplicates_recalibrate {
	tag "${meta.sampleName}"

	input:
		val (reference)
		val (common_vcf)
		tuple val (meta), val (type), path (bam)

	script:
	"""#!/usr/bin/env bash
mkdir temp
mkdir QC
# next give it to SAMBLASTER for read dup marking, followed by sorting back to coordinates and format to BAM
/opt/bin/sambamba sort \\
	--sort-by-name \\
	-t ${params.sambamba.threads} \\
	-m ${params.sambamba.ram}GB \\
	--tmpdir temp \\
	-o /dev/stdout \\
	${bam} | samtools view -h | /opt/samblaster-0.1.26/samblaster | samtools view -Sb | /opt/bin/sambamba sort \\
	-t ${params.sambamba.threads} \\
	-m ${params.sambamba.ram}GB \\
	--tmpdir temp \\
	-o ${meta.sampleName}.${type}.cleaned.sorted.marked.bam \\
	/dev/stdin

java -Xmx${params.picard.ram}G -Dpicard.useLegacyParser=false \\
		-jar ${params.picard.jar} AddOrReplaceReadGroups \\
		-I ${meta.sampleName}.${type}.cleaned.sorted.marked.bam \\
		-O ${meta.sampleName}.${type}.cleaned.sorted.readgroups.marked.bam \\
		-ID 1 -LB Lib1 -PL ILLUMINA -PU Run1 -SM $type \\
		-MAX_RECORDS_IN_RAM ${params.picard.max_records_in_ram}

java -Xmx${params.gatk.ram}G -jar ${params.gatk.jar} BaseRecalibrator \\
	-R ${reference} \\
	-I ${meta.sampleName}.${type}.cleaned.sorted.readgroups.marked.bam \\
	--known-sites ${common_vcf} \\
	--use-original-qualities \\
	-O QC/${meta.sampleName}.${type}.GATK4.pre.recal.table

java -Xmx${params.gatk.ram}G -jar ${params.gatk.jar} ApplyBQSR \\
	-R ${reference} \\
	-I ${meta.sampleName}.${type}.cleaned.sorted.readgroups.marked.bam \\
	-O ${meta.sampleName}.${type}.bam \\
	-bqsr QC/${meta.sampleName}.${type}.GATK4.pre.recal.table

	rm ${meta.sampleName}.${type}.cleaned.sorted.readgroups.marked.bam
	rm ${meta.sampleName}.${type}.cleaned.sorted.readgroups.marked.bam.bai

java -Xmx${params.gatk.ram}G -jar ${params.gatk.jar} BaseRecalibrator \\
	-R ${reference} \\
	-I ${meta.sampleName}.${type}.bam \\
	--known-sites ${common_vcf} \\
	--use-original-qualities \\
	-O QC/${meta.sampleName}.${type}.GATK4.post.recal.table

/opt/bin/sambamba index \\
	-t ${params.sambamba.threads} \\
	${meta.sampleName}.${type}.bam \\
	${meta.sampleName}.${type}.bam.bai

	"""

}


