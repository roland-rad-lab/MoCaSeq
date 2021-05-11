




process cnv_kit_matched {
	tag "${meta.sampleName}"


	input:
		val (reference)
		val (interval_bed)
		val (meta), path (bam_normal), path (bai_normal), path (bam_tumor), path (bai_tumor)

	script:
	"""#!/usr/bin/env bash

	echo "Hello from CNVKit"

	"""

}


