process delly_matched_call {
	 tag "${meta.sampleName}"

	input:
	val (reference)
	tuple val(meta), path (bam_normal), path (bam_tumor)

	output:
	tuple val(meta), path("${meta.sampleName}.pre.bcf"), emit: result

	script:
	"""#!/usr/bin/env bash
delly call \
	-o ${meta.sampleName}.pre.bcf \
	-g ${reference} \
	${bam_tumor} \
	${bam_normal}
	"""
}

process delly_matched_filter {
	 tag "${meta.sampleName}"

	input:
	val (reference)
	tuple val(meta), path (delly_pre_bcf)

	output:
	tuple val(meta), path("${meta.sampleName}.delly.vcf"), emit: result

	script:
	"""#!/usr/bin/env bash
delly filter \
	-f somatic \
	-o ${meta.sampleName}.delly.bcf \
	-s Samples.tsv \
	${delly_pre_bcf}

bcftools view ${meta.sampleName}.delly.bcf > ${meta.sampleName}.delly.vcf
	"""
}

