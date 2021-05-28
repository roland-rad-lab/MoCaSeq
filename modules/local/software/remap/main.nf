
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
		path ("QC/*")

	script:
	"""#!/usr/bin/env bash

fastqc -t ${params.fastqc.threads} \\
	${fastq_r1} \\
	${fastq_r2} \\
	--outdir=QC

	"""
}

