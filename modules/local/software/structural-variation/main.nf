



process structural_variation_matched {
	tag "${meta.sampleName}"

	input:
		tuple path (interval_bed), path(interval_bed_index)
		tuple val (meta)


	script:
	"""#!/usr/bin/env bash


	"""
}



