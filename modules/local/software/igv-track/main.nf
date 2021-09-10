
process igv_track_depth {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Tracks", mode: "copy"

	input:
		val (genome_build)
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

process igv_track_cnr {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Tracks", mode: "copy"

	input:
		val (genome_build)
		tuple path (interval_bed), path (interval_bed_index)
		tuple val (meta), val (type), val(coverage_source), path (cnr)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.${coverage_source}.bigWig")

	script:
	"""#!/usr/bin/env Rscript

library (dplyr)

data <- read.table (file="${cnr}",sep="\\t",header=T,stringsAsFactors=F)
data_genome_sizes <- read.table (file=gzfile ("${interval_bed}"),sep="\\t",header=F,stringsAsFactors=F)
names (data_genome_sizes) <- c("chromosome", "start", "end")

head (data)
head (data_genome_sizes)

write.table (data_genome_sizes %>% dplyr::select (chromosome,end) %>% data.frame ,file="genome.sizes",sep="\\t",row.names=F,col.names=F,quote=F)

data_bed <- switch ("${coverage_source}",
		"dryclean"={
			data %>%
			dplyr::filter (chromosome %in% !!data_genome_sizes[,"chromosome"]) %>%
			dplyr::select (chromosome,start,end,log2)
		},
		{
			data
		}) %>%
		dplyr::mutate (start=as.numeric (start)) %>%
		dplyr::arrange (chromosome,start) %>%
		data.frame

write.table (data_bed,file="${meta.sampleName}.${type}.${coverage_source}.bedGraph",sep="\\t",row.names=F,col.names=F,quote=F)

system ("bedGraphToBigWig ${meta.sampleName}.${type}.${coverage_source}.bedGraph genome.sizes ${meta.sampleName}.${type}.${coverage_source}.bigWig")
file.remove ("${meta.sampleName}.${type}.${coverage_source}.bedGraph")

	"""
}

process igv_track_cns {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Tracks", mode: "copy"

	input:
		val (genome_build)
		tuple val (meta), val (type), val (coverage_source), path (cns)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.${coverage_source}.bedGraph")

	script:
	"""#!/usr/bin/env Rscript

library (dplyr)

data <- read.table (file="${cns}",sep="\\t",header=T,stringsAsFactors=F)
#head (data)

data_bed <- data %>%
	dplyr::select (chromosome,start,end,log2) %>%
	data.frame

write.table (data_bed,file="${meta.sampleName}.${type}.${coverage_source}.bedGraph",sep="\\t",quote=F,row.names=F)

	"""
}

process igv_track_rds {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Tracks", mode: "copy"

	input:
		val (genome_build)
		tuple path (interval_bed), path (interval_bed_index)
		val (coverage_source)
		tuple val (meta), val (type), val (resolution), path (coverage_rds)

	output:
		tuple val (meta), val (type), val (resolution), path ("${meta.sampleName}.${type}.${coverage_source}.${resolution}.bigWig")

	script:
	"""#!/usr/bin/env Rscript
library (dplyr)

interval_file <- gzfile ("${interval_bed}", 'rt')
data_interval <- read.table (file=interval_file,sep="\\t",header=F,stringsAsFactors=F)
names (data_interval) <- c("Chrom", "Start", "End")
head (data_interval)

write.table (data_interval %>% dplyr::mutate (Size=End-Start) %>% dplyr::select (Chrom,Size) %>% data.frame,file="intervals.sizes",sep="\\t",row.names=F,col.names=F,quote=F)

data <- as.data.frame (readRDS ("${coverage_rds}")) %>%
	dplyr::filter (!is.na(reads.corrected)) %>%
	dplyr::select (seqnames,start,end,reads.corrected) %>%
	dplyr::mutate (seqnames=as.character(seqnames)) %>%
	dplyr::arrange (seqnames,start) %>%
	data.frame

write.table (data,file="${meta.sampleName}.${type}.${coverage_source}.${resolution}.bedGraph",row.names=F,col.names=F,sep="\\t",quote=F)
system ("bedGraphToBigWig ${meta.sampleName}.${type}.${coverage_source}.${resolution}.bedGraph intervals.sizes ${meta.sampleName}.${type}.${coverage_source}.${resolution}.bigWig")
system ("rm ${meta.sampleName}.${type}.${coverage_source}.${resolution}.bedGraph")
	"""
}


