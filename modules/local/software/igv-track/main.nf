
process igv_track_depth {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${meta.sampleName}/results/Tracks", mode: "copy", pattern: "*.bigWig"

	input:
		tuple val (meta), path (bam), path (bai)
		val (intervals)
		path (interval_bed)

	output:
		tuple val (meta), path ("${meta.sampleName}.depth.raw.all.bigWig")

	script:
	"""#!/usr/bin/env bash
source ${params.script_base}/file_handling.sh
temp_file=\$(moc_mktemp_file .)
trap "rm \${temp_file}" EXIT

extract_if_zip ${interval_bed} interval_bed_extracted \${temp_file}

cat \${interval_bed_extracted} | awk 'BEGIN{FS=OFS="\\t";}{print \$1,\$3-\$2;}' > genome.sizes

for chromosome in ${intervals};
do
	echo "fixedStep chrom=\${chromosome} start=1 step=1" >> ${meta.sampleName}.depth.raw.all.wig
	samtools depth -ar \${chromosome} ${bam} | cut -f 3) >> ${meta.sampleName}.depth.raw.all.wig
done
wigToBigWig ${meta.sampleName}.depth.raw.all.wig genome.sizes ${meta.sampleName}.depth.raw.all.bigWig
rm ${meta.sampleName}.depth.raw.all.wig

	"""
}


