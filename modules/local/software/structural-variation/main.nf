



process structural_variation_matched {
	tag "${meta.sampleName}"

	input:
		tuple path (interval_bed), path(interval_bed_index)
		tuple val (meta), path (manta_vcf), path (manta_vcf_index), path (delly_bcf), path (normal_cns), path (tumor_cns)


	script:
	"""#!/usr/bin/env bash


	"""
}



