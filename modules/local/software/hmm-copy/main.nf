
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
		tuple val (resolution), val (meta), path (normal_wig), path (tumor_wig), val (gc_wig), val (map_wig)



	script:
	"""#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(HMMcopy))
suppressPackageStartupMessages(library(DNAcopy))
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(naturalsort))
suppressPackageStartupMessages(library(GenomicRanges))

# read in wig files and correct for GC and mappability bias
normal <- wigsToRangedData("${normal_wig}","${gc_wig}","${map_wig}")
normal\$reads <- normal\$reads+1
normal <- as.data.frame(correctReadcount(normal))
normal_copy=GRanges(normal\$chr, IRanges(normal\$start, normal\$end),copy=normal\$copy)

tumor <- wigsToRangedData("${tumor_wig}","${gc_wig}","${map_wig}")
tumor\$reads <- tumor\$reads+1
tumor <- as.data.frame(correctReadcount(tumor))
tumor_copy=GRanges(tumor\$chr, IRanges(tumor\$start, tumor\$end),copy=tumor\$copy)

	"""
}



