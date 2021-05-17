
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

	output:
		tuple val (meta), val (resolution), path ("${meta.sampleName}.HMMCopy.${resolution}.log2RR.txt"), path ("${meta.sampleName}.HMMCopy.${resolution}.segments.txt"), emit: result

	script:

		def flank_length = params.hmm_copy.flank_length["${meta.organism}"]
		def filter_regions = params.hmm_copy.filter_regions["${meta.organism}"]

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

  cat(paste("track type=wiggle_0 name=\\"", sample, "\\"", sep = ""),
    file = file, sep = "\\n")
  temp <- data.frame(chr = correctOutput\$chr, dat)
  width <- correctOutput\$start[2] - correctOutput\$start[1]
  chrs <- levels(correctOutput\$chr)

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
      " span=", width, sep = ""), file = file, append = TRUE, sep = "\\n")
    cat(out, file = file, sep = "\\n", append = TRUE)
  }
}
# ^Had to rewrite that becuase subset(temp, chr == chr) is nonsense

intervals <- strsplit ("${intervals}", ",", fixed=T)[[1]]

gc_ranged <- wigToRangedData ("${gc_wig}")
setkeyv(gc_ranged,"chr")
gc_ranged_subset <- gc_ranged[.(intervals)]

gc_chr_level_info <- gc_ranged_subset[,.(count=.N),by=chr]
gc_chr_levels <- gc_chr_level_info[order(gc_chr_level_info[,"count"],decreasing=T),"chr"]

rangedDataToWig (gc_ranged_subset[,chr:=factor (chr,levels=gc_chr_levels\$chr,ordered=T)], "intervals.gc.wig",column="value")
rm (gc_ranged_subset, gc_ranged, gc_chr_levels, gc_chr_level_info)

map_ranged <- wigToRangedData ("${map_wig}")
setkeyv(map_ranged,"chr")
map_ranged_subset <- map_ranged[.(intervals)]

map_chr_level_info <- map_ranged_subset[,.(count=.N),by=chr]
map_chr_levels <- map_chr_level_info[order(map_chr_level_info[,"count"],decreasing=T),"chr"]

rangedDataToWig (map_ranged_subset[,chr:=factor (chr,levels=map_chr_levels\$chr,ordered=T)], "intervals.map.wig", column="value")
rm (map_ranged_subset, map_ranged, map_chr_levels, map_chr_level_info)


# read in wig files and correct for GC and mappability bias
normal <- wigsToRangedData("${normal_wig}","intervals.gc.wig","intervals.map.wig")
normal\$reads <- normal\$reads+1
normal <- as.data.frame(correctReadcount(normal))
normal_copy=GRanges(normal\$chr, IRanges(normal\$start, normal\$end),copy=normal\$copy)

tumor <- wigsToRangedData("${tumor_wig}","intervals.gc.wig","intervals.map.wig")
tumor\$reads <- tumor\$reads+1
tumor <- as.data.frame(correctReadcount(tumor))
tumor_copy=GRanges(tumor\$chr, IRanges(tumor\$start, tumor\$end),copy=tumor\$copy)

flankLength=${flank_length}
filterRegions=read.delim ("${filter_regions}")

colnames(filterRegions)[1:3] <- c("space","start","end")
filterRegions\$start <- filterRegions\$start - flankLength
filterRegions\$end <- filterRegions\$end + flankLength
filtering=GRanges(filterRegions\$space, IRanges(filterRegions\$start, filterRegions\$end))

hits <- findOverlaps(query = normal_copy, subject = filtering)
ind <- queryHits(hits)
message("Removed ", length(ind), " bins near centromeres.")
normal_copy=(normal_copy[-ind, ])

hits <- findOverlaps(query = tumor_copy, subject = filtering)
ind <- queryHits(hits)
message("Removed ", length(ind), " bins near centromeres.")
tumor_copy=(tumor_copy[-ind, ])

# computation of the copy number states from the log fold change
somatic_copy <- tumor_copy
somatic_copy\$copy <- tumor_copy\$copy - normal_copy\$copy
somatic_tab <- as.data.frame(somatic_copy)
colnames(somatic_tab) <- c("Chrom", "Start", "End", "width", "strand", "log2Ratio")
somatic_tab <- somatic_tab[,c("Chrom", "Start", "End","log2Ratio")]
write.table(somatic_tab,"${meta.sampleName}.HMMCopy.${resolution}.log2RR.txt",quote=F,row.names=F,col.names=T,sep="\\t")

# segmentation of the CN plot
somatic_CNA <- smooth.CNA(CNA(genomdat=somatic_tab\$log2Ratio,chrom=somatic_tab\$Chrom,maploc=somatic_tab\$Start,data.type="logratio"))

cnv_segments <- segment(somatic_CNA,alpha=0.0001,min.width=5,undo.splits='sdundo',undo.SD=2,verbose=2)\$output

colnames(cnv_segments) <- c("ID", "Chrom", "Start", "End", "num.mark", "Mean")
cnv_segments <- cnv_segments[,c("Chrom", "Start", "End","Mean")]
cnv_segments <- cnv_segments[naturalorder(cnv_segments\$Chrom),]
write.table(cnv_segments,"${meta.sampleName}.HMMCopy.${resolution}.segments.txt",quote=F,row.names=F,col.names=T,sep="\\t")

	"""
}

process hmm_copy_plot {
	tag "${meta.sampleName}"

	input:
		tuple path (interval_bed), path(interval_bed_index)
		tuple val (meta), val (resolution), path (log2_file), path (segments_file)

	script:
	"""#!/usr/bin/env Rscript
library (dplyr)

interval_file <- gzfile ("${interval_bed}", 'rt')
data_interval <- read.table (file=interval_file,sep="\\t",header=F,stringsAsFactors=F)
names (data_interval) <- c("Chrom", "Start", "End")
head (data_interval)

data_ratio <- read.table (file="${log2_file}",sep="\\t",header=T,stringsAsFactors=F)
head (data_ratio)

data_segments <- read.table (file="${segments_file}",sep="\\t",header=T,stringsAsFactors=F)
head (data_segments)

data_interval_plot <- data_interval %>%
	mutate (CumulativeStart=cumsum (End)-End) %>%
	mutate (CumulativeMidpoint=CumulativeStart+((Start+End)/2)) %>%
	data.frame

 data_ratio_plot <- data_ratio %>%
filter (!is.na (log2Ratio)) %>%
inner_join (data_interval_plot,by="Chrom",suffix=c("",".Chrom")) %>%
mutate (Start.Genome=Start+CumulativeStart,End.Genome=End+CumulativeStart) %>%
data.frame

data_segments_plot <- data_segments %>%
inner_join (data_interval_plot,by="Chrom",suffix=c("",".Chrom")) %>%
mutate (Start.Genome=Start+CumulativeStart,End.Genome=End+CumulativeStart) %>%
data.frame

head (data_interval_plot)
head (data_ratio_plot)
head (data_segments_plot)


pdf (file="genome.pdf",width=9,height=6)

ggplot (data_ratio_plot) +
	geom_segment (aes(x=Start.Genome,y=log2Ratio,xend=End.Genome,yend=log2Ratio)) +
	geom_segment (data=data_segments_plot,aes(x=Start.Genome,y=Mean,xend=End.Genome,yend=Mean),colour="red") +
	geom_vline (data=data_interval_plot,aes(xintercept=CumulativeStart)) +
	geom_text (data=data_interval_plot,aes(x=CumulativeMidpoint,y=2,label=Chrom)) +
	ylim (-2,2) +
	theme_bw () +
	theme (
		panel.grid.major.x=element_blank (),
		panel.grid.minor.x=element_blank (),
		axis.ticks.x=element_blank (),
		axis.text.x=element_blank ()
	)

dev.off ()

chromosomes <- data_interval %>% pull (Chrom)
plot_list <- vector ("list",length (chromosomes))


for ( i in seq_along (chromosomes) )
{
	plot_list[[i]] <- ggplot (data_ratio_plot %>% filter (Chrom==!!chromosomes[[i]]) %>% data.frame) +
		geom_segment (aes(x=Start,y=log2Ratio,xend=End,yend=log2Ratio)) +
		geom_segment (data=data_segments_plot %>% filter (Chrom==!!chromosomes[[i]]) %>% data.frame,aes(x=Start,y=Mean,xend=End,yend=Mean),colour="red") +
		ylim (-2,2) +
		labs (title=chromosomes[[i]]) +
		theme_bw () +
		theme (
			panel.grid.major.x=element_blank (),
			panel.grid.minor.x=element_blank (),
			#axis.ticks.x=element_blank (),
			axis.text.x=element_blank (),
			#plot.margin = unit(c(4,1,1,1), "cm")
		)
}

# Need to do this outside of pdf call to prevent blank first page
p <- marrangeGrob (plot_list,nrow=1,ncol=1) 

pdf (file="chromosomes.pdf",width=9)
p
dev.off ()


	"""
}


