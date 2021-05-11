




process cnv_kit_matched {
	tag "${meta.sampleName}"


	input:
		val (reference)
		val (reference_flat)
		tuple val (interval_bed), val (interval_bed_index)
		tuple val (meta), path (bam_normal), path (bai_normal), path (bam_tumor), path (bai_tumor)

	script:
	"""#!/usr/bin/env bash
cnvkit.py batch \\
	${bam_tumor} \\
	--normal ${bam_normal} \\
	--fasta ${reference} \\
	--output-reference Reference.cnn \\
	--output-dir . \\
	--short-names \\
	--diagram \\
	--scatter \\
	--annotate ${reference_flat} \\
	--access "" \\
	--targets ${interval_bed} \\
	--drop-low-coverage \\
	-m wgs \\
	-p ${params.cnv_kit.threads}"

	"""

}


