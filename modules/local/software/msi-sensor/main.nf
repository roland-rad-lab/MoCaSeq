

params.msi = [:]

process msi_matched {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/msisensor", mode: "copy"

	input:
		val (genome_build)
		val (micro_satellite)
		tuple val (meta), path (bam_normal), path (bai_normal), path (bam_tumor), path (bai_tumor)

	output:
		tuple val (meta), path ("${meta.sampleName}.msisensor"), emit: result

	script:
	"""#!/usr/bin/env bash

msisensor msi \\
	-n ${bam_normal} \\
	-t ${bam_tumor} \\
	-o ${meta.sampleName}.msisensor \\
	-d ${micro_satellite} \\
	-b ${params.msi.threads}
	"""

	stub:
	"""#!/usr/bin/env bash
#cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/msisensor/${meta.sampleName}.msisensor .
touch ${meta.sampleName}.msisensor
	"""
}


