process frag_counter_wig_to_rds {

	publishDir "${params.output_base}/PON", mode: "copy"

	input:
		val (genome_name)
		val (intervals)
		tuple val (resolution), val (gc_wig), val (map_wig)

	output:
		tuple val (resolution), path ("${genome_name}.gc.${resolution}.rds"), path ("${genome_name}.map.${resolution}.rds"), emit: result

	script:
	"""#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(HMMcopy))
suppressPackageStartupMessages(library(GenomicRanges))

intervals <- strsplit ("${intervals}", ",", fixed=T)[[1]]

gc_ranged <- wigToRangedData ("${gc_wig}")
setkeyv(gc_ranged,"chr")
gc_ranged_subset <- gc_ranged[.(intervals)]

gc_chr_level_info <- gc_ranged_subset[,.(count=.N),by=chr]
gc_chr_levels <- gc_chr_level_info[order(gc_chr_level_info[,"count"],decreasing=T),"chr"]
gc_ranged_subset[,c("chr", "end") := list (factor (chr,levels=gc_chr_levels\$chr,ordered=T), end -1 )]
data.table::setnames (gc_ranged_subset, "value", "score")

saveRDS (GenomicRanges::makeGRangesFromDataFrame (as.data.frame (gc_ranged_subset),keep.extra.columns=T),file="${genome_name}.gc.${resolution}.rds")
rm (gc_ranged_subset, gc_ranged, gc_chr_levels, gc_chr_level_info)

map_ranged <- wigToRangedData ("${map_wig}")
setkeyv(map_ranged,"chr")
map_ranged_subset <- map_ranged[.(intervals)]

map_chr_level_info <- map_ranged_subset[,.(count=.N),by=chr]
map_chr_levels <- map_chr_level_info[order(map_chr_level_info[,"count"],decreasing=T),"chr"]
map_ranged_subset[,c("chr", "end") := list (factor (chr,levels=map_chr_levels\$chr,ordered=T), end -1)]
data.table::setnames (map_ranged_subset, "value", "score")

saveRDS (GenomicRanges::makeGRangesFromDataFrame (as.data.frame (map_ranged_subset),keep.extra.columns=T),file="${genome_name}.map.${resolution}.rds")
	"""

	stub:
	"""#!/usr/bin/env bash
#cp ${params.stub_dir}/PON/${genome_name}.gc.${resolution}.rds .
#cp ${params.stub_dir}/PON/${genome_name}.map.${resolution}.rds .
touch ${genome_name}.gc.${resolution}.rds
touch ${genome_name}.map.${resolution}.rds
	"""
}

process frag_counter {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${meta.sampleName}/results/fragCounter", mode: "copy"

	input:
		tuple val (resolution), path ("gc${resolution}.rds"), path ("map${resolution}.rds")
		tuple val (meta), val (type), path (bam), path (bai)

	output:
		tuple val (meta), val (type), val (resolution), path ("${meta.sampleName}.${type}.coverage.${resolution}.rds"), emit: result

	script:
	"""#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(fragCounter))
suppressPackageStartupMessages(library(skidb))
suppressPackageStartupMessages(library(data.table))

out = fragCounter(bam="${bam}", window=${resolution}, gc.rds.dir=".", map.rds.dir=".")
saveRDS (out,file="${meta.sampleName}.${type}.coverage.${resolution}.rds")
	"""

	stub:
	"""#!/usr/bin/env bash
cp ${params.stub_dir}/${meta.sampleName}/results/fragCounter/${meta.sampleName}.${type}.coverage.${resolution}.rds .
#touch ${meta.sampleName}.${type}.coverage.${resolution}.rds
	"""
}

