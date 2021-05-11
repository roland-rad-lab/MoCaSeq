
params.manta = [:]

process manta_matched {
	tag "${meta.sampleName}"

	input:
		val (reference)
		val (interval_bed)
		tuple val (meta), path (bam_normal), path (bai_normal), path (bam_tumor), path (bai_tumor)

	output:
		tuple val (meta), val ("matched"), path ("results/variants/somaticSV.vcf.gz"), path ("results/variants/somaticSV.vcf.gz.tbi"), emit: sv
		tuple val (meta), val ("matched"), path ("results/variants/candidateSmallIndels.vcf.gz"), path ("results/variants/candidateSmallIndels.vcf.gz.tbi"), emit: indel

	script:
	"""#!/usr/bin/env bash
tabix -p bed ${interval_bed}
python2 ${params.manta.dir}/bin/configManta.py \\
	--normalBam ${bam_normal} \\
	--tumorBam ${bam_tumor} \\
	--referenceFasta ${reference} \\
	--runDir . \\
	--callRegions ${interval_bed} \\
	--generateEvidenceBam

python2 runWorkflow.py -m local -j ${params.manta.threads}

	"""
}

process manta_matched_post {
	tag "${meta.sampleName}"

	input:
		tuple val (meta), val(type), path (somatic_vcf), path (somatic_vcf_index)

	output:
		tuple val (meta), path ("${meta.sampleName}.Manta.annotated.vcf"), emit: result

	script:

	def snpeff_version = params.annotation.snpeff["${meta.organism}"]

	"""#!/usr/bin/env bash

zcat ${somatic_vcf} | java -jar ${params.annotation.snpeff.dir}/SnpSift.jar filter "( ( ( FILTER = 'PASS' ) & ( ( GEN[Tumor].SR[1] + GEN[Tumor].SR[0] ) * 0.05 <= GEN[Tumor].SR[1] ) ) | ( ( FILTER = 'PASS' ) & ( ( GEN[Tumor].PR[1] + GEN[Tumor].PR[0] ) * 0.05 <= GEN[Tumor].PR[1] ) ) )" | bgzip -c > ${meta.sampleName}.Manta.vcf.gz
tabix -p vcf ${meta.sampleName}.Manta.vcf.gz
bcftools stats ${meta.sampleName}.Manta.vcf.gz > ${meta.sampleName}.Manta.vcf.gz.stats

java -Xmx${params.annotation.snpeff.ram}g \\
	-jar ${params.annotation.snpeff.dir}/snpEff.jar \\
	${snpeff_version} -canon \\
	-csvStats ${meta.sampleName}.Manta.annotated.vcf.stats \\
	${meta.sampleName}.Manta.vcf.gz \\
	> ${meta.sampleName}.Manta.annotated.vcf

	"""
}





