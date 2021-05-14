
process hmm_copy_wig {
	tag "${meta.sampleName}"


	input:
		val (intervals)
		each (resolution)
		tuple val (meta), val (type), path (bam), path (bai)

	output:
		tuple val (meta), val (type), val (resolution), path ("HMMCopy/${meta.sampleName}.${type}.${resolution}.wig"), emit: result

	script:
	"""#!/usr/bin/env bash

mkdir HMMCopy
${params.hmm_copy.dir}/bin/readCounter -w ${resolution} -q20 -c ${intervals} ${bam} > HMMCopy/${meta.sampleName}.${type}.${resolution}.wig

	"""
}

