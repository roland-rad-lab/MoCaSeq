
process igv_track_depth {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${meta.sampleName}/results/Tracks", mode: "copy", pattern: "*.bigWig"

	input:
		val (intervals)
		tuple path (interval_bed), path (interval_bed_index)
		tuple val (meta), val (type), path (bam), path (bai)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.depth.raw.all.bigWig")

	script:
	"""#!/usr/bin/env bash

zcat ${interval_bed} | awk 'BEGIN{FS=OFS="\\t";}{print \$1,\$3-\$2;}' > genome.sizes

for chromosome in ${intervals};
do
	echo "fixedStep chrom=\${chromosome} start=1 step=1" >> ${meta.sampleName}.${type}.depth.raw.all.wig
	samtools depth -ar \${chromosome} ${bam} | cut -f 3 >> ${meta.sampleName}.${type}.depth.raw.all.wig
done
wigToBigWig ${meta.sampleName}.${type}.depth.raw.all.wig genome.sizes ${meta.sampleName}.${type}.depth.raw.all.bigWig
rm ${meta.sampleName}.${type}.depth.raw.all.wig

	"""
}


