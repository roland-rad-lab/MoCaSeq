

process mutect_matched {
	tag "${meta.sampleName}"

	input:
		tuple val (meta), path (loh_tsv), path (segments_tsv)


	script:
	"""#!/usr/bin/env Rscript

library(ggplot2)
library (dplyr)
library(BubbleTree)


# function to plot all chromosomes individually
PlotChromBT <- function(data=pred, chrom, mode="adjusted"){
  
  # define the data, chromosome and output name for both modes
  if(mode=="adjusted"){
    dat <- pred@rbd.adj
    plotdata <- dat[dat$seqnames == as.character(chrom),]
    pdfname <- paste0(name,"/results/BubbleTree/",name,"_Chromosomes/adjusted/",name,".Bubbletree.adjusted.chr",chrom,".pdf")
  } else if(mode == "unadjusted"){
    dat <- pred@rbd
    plotdata <- dat[dat$seqnames == as.character(chrom),]
    pdfname <- paste0(name,"/results/BubbleTree/",name,"_Chromosomes/unadjusted/",name,".Bubbletree.unadjusted.chr",chrom,".pdf")
  } else {
    stop("ERROR: only adjusted and unadjusted valid for mode")
  }

  # for some stupid reason, the first entry is set to color white, this way we catch it
  plotdata <- rbind(emptyElem, plotdata)

  # draw the plot
  btree <- drawBTree(btreeplotter, plotdata)
  
  # save the plot
  ggsave(filename = pdfname, plot = btree, width = 16, height = 9, device = "pdf")
}

data_loh <- read.table (file="${loh_tsv}",sep="\\t",header=T,stringsAsFactors=F)
data_cnv <- read.table (file="${segments_tsv}",sep="\\t",header=T,stringsAsFactors=F)

gr_loh <- data_loh %>%
	dplyr::mutate (start=Pos) %>%
	dplyr::rename () %>%
	GenomicRanges ()

gr_cnv <- data_cnv %>%
	dplyr::rename () %>%
	GenomicRanges ()

print (gr_loh)
print (gr_cnv)



# for some stupid reason, the first entry is set to color white, this way we catch it
emptyElem <- data.frame(seqnames=0,start=0, end=0,width=0, strand="*", seg.id=1000, num.mark=1000, lrr=0, kurtosis=0, hds=0, hds.sd=0,het.cnt=0,seg.size=0)

	"""

}



