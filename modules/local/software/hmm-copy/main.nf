
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

process hmm_copy_tsv {
	tag "${meta.sampleName}"

	input:
		val (intervals)
		tuple val (resolution), val (meta), path (normal_wig), path (tumor_wig), val (gc_wig), val (map_wig)

	script:
	"""#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(HMMcopy))
suppressPackageStartupMessages(library(DNAcopy))
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(naturalsort))
suppressPackageStartupMessages(library(GenomicRanges))

rangedDataToWig <- function(correctOutput, file, column = "copy", sample = "R", verbose = TRUE) {
  dat <- c(correctOutput[[column]])
  if (length(dat) == 0) {
    stop(paste(column, "is not a valid column"))
  }
  dat[is.na(dat)] <- -1

  cat(paste("track type=wiggle_0 name=\"", sample, "\"", sep = ""),
    file = file, sep = "\n")
  temp <- data.frame(chr = correctOutput$chr, dat)
  width <- correctOutput$start[2] - correctOutput$start[1]
  chrs <- levels(correctOutput$chr)

  for (i in 1:length(chrs)) {
    #chr <- chrs[i]
    #out <- subset(temp, chr == chr)[, 2]
    not_that_chr <- chrs[i]
    out <- subset (temp, chr == not_that_chr)[,2]
    if (verbose) {
      message(paste("Outputting chromosome ", not_that_chr,
        " (", length(out), ")", sep = ""))
    }
    cat(paste("fixedStep chrom=", not_that_chr, " start=1 step=", width,
      " span=", width, sep = ""), file = file, append = TRUE, sep = "\n")
    cat(out, file = file, sep = "\n", append = TRUE)
  }
}
# ^Had to rewrite that becuase subset(temp, chr == chr) is nonsense

intervals <- strsplit ("${intervals}", ",", fixed=T)[[1]]

gc_ranged <- wigToRangedData ("${gc_wig}")
setkeyv(gc_ranged,"chr")
gc_ranged_subset <- gc_ranged[.(intervals)]

gc_chr_level_info <- gc_ranged_subset[,.(count=.N),by=chr]
gc_chr_levels <- gc_chr_level_info[order(gc_chr_level_info[,"count"],decreasing=T),"chr"]

rangedDataToWig (gc_ranged_subset[,chr:=factor (chr,levels=gc_chr_levels$chr,ordered=T)], "intervals.gc.wig",column="value")
rm (gc_ranged_subset, gc_ranged, gc_chr_levels, gc_chr_level_info)

map_ranged <- wigToRangedData ("${map_wig}")
setkeyv(map_ranged,"chr")
map_ranged_subset <- map_ranged[.(intervals)]

map_chr_level_info <- map_ranged_subset[,.(count=.N,by=chr]
map_chr_levels <- map_chr_level_info[order(map_chr_level_info[,"count"],decreasing=T),"chr"]

rangedDataToWig (map_ranged_subset[,chr:=factor (chr,levels=map_chr_levels$chr,ordered=T)], "intervals.map.wig", column="value")
rm (map_ranged_subset, map_ranged, map_chr_levels, map_chr_level_info)


# read in wig files and correct for GC and mappability bias
normal <- wigsToRangedData("${normal_wig}","intervals.gc.wig","intervals.map.wig")
normal\$reads <- normal\$reads+1
normal <- as.data.frame(correctReadcount(normal))
normal_copy=GRanges(normal\$chr, IRanges(normal\$start, normal\$end),copy=normal\$copy)

tumor <- wigsToRangedData("${tumor_wig}","intevals.gc.wig","intervals.map.wig")
tumor\$reads <- tumor\$reads+1
tumor <- as.data.frame(correctReadcount(tumor))
tumor_copy=GRanges(tumor\$chr, IRanges(tumor\$start, tumor\$end),copy=tumor\$copy)

	"""
}



