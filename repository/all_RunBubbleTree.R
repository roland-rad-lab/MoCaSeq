#!/usr/bin/Rscript

##########################################################################################
##
## all_RunBubbleTree.R
##
## Run BubbleTree.
##
##########################################################################################

#parallel 'Rscript /media/rad/SSD1/DNA/repository/all_RunBubbleTree.R {} Copywriter' ::: $(ls)

args <- commandArgs(TRUE)

name = args[1]
cnv_method = args[2]

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(BubbleTree))
suppressPackageStartupMessages(library(GenomicRanges))

# for some stupid reason, the first entry is set to color white, this way we catch it
emptyElem <- data.frame(seqnames=0,start=0, end=0,width=0, strand="*", seg.id=1000, num.mark=1000, lrr=0, kurtosis=0, hds=0, hds.sd=0,het.cnt=0,seg.size=0)

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

loh = read.delim(paste0(name,"/results/LOH/",name,".VariantsForLOH.txt"))
snp.gr <- GenomicRanges::GRanges(loh$Chrom, IRanges(loh$Pos, loh$Pos), freq=loh[,"Tumor_Freq"])

if (cnv_method=="Copywriter"){
	cnv = read.delim(paste0(name,"/results/Copywriter/",name,".seg.dat.fn"))
} else if (cnv_method=="HMMCopy") {
	cnv = read.delim(paste0(name,"/results/HMMCopy/",name,".HMMCopy.20000.segments.txt.fn"))
}

cnv.gr <- GenomicRanges::GRanges(cnv$Chromosome, IRanges(cnv$Start, cnv$End), num.mark=cnv$Num_Probes, seg.mean=cnv$Segment_Mean)

r <- new("RBD", unimodal.kurtosis=-0.1)

rbd <- makeRBD(r, snp.gr,cnv.gr)
btreepredictor <- new("BTreePredictor", rbd=rbd, max.ploidy=6, prev.grid=seq(0.2,3, by=0.01))
btreepredictor@config$min.segSize <- ifelse(max(btreepredictor@rbd$seg.size,na.rm=TRUE) < 0.4, 0.1, 0.4)
pred <- btpredict(btreepredictor)


# plot full unadjusted
btreeplotter <- new("BTreePlotter", branch.col="gray50")
pred@rbd <- rbind(emptyElem, pred@rbd) # for some stupid reason, the first entry is set to color white, this way we catch it
btree <- drawBTree(btreeplotter, pred@rbd)
pdfname <- paste0(name,"/results/BubbleTree/",name,".Bubbletree.unadjusted.pdf")
ggsave(filename = pdfname, plot = btree, width = 16, height = 9, device = "pdf")

# plot full adjusted
btreeplotter <- new("BTreePlotter", branch.col="gray50")
pred@rbd.adj <- rbind(emptyElem, pred@rbd.adj) # for some stupid reason, the first entry is set to color white, this way we catch it
btree <- drawBTree(btreeplotter, pred@rbd.adj)
pdfname <- paste0(name,"/results/BubbleTree/",name,".Bubbletree.adjusted.pdf")
ggsave(filename = pdfname, plot = btree, width = 16, height = 9, device = "pdf")

# plot shift from unadjusted to adjusted
btreeplotter <- new("BTreePlotter", branch.col="gray50")
rbd1 <- pred@rbd
rbd2 <- pred@rbd.adj
arrows <- trackBTree(btreeplotter, rbd1, rbd2, min.srcSize=0.01, min.trtSize=0.01)
btree <- drawBTree(btreeplotter, rbd2) + arrows
pdfname <- paste0(name,"/results/BubbleTree/",name,".Bubbletree.adjustmentshift.pdf")
ggsave(filename = pdfname, plot = btree, width = 16, height = 9, device = "pdf")


# plot all single chromosomes for both

# create output dirs
dir.create(paste0(name,"/results/BubbleTree/",name,"_Chromosomes/adjusted"), showWarnings = FALSE, recursive = T)
dir.create(paste0(name,"/results/BubbleTree/",name,"_Chromosomes/unadjusted"), showWarnings = FALSE, recursive = T)

# generate a plot for each chromosome
allchroms <- as.character(unique(pred@rbd$seqnames))
allchroms <- allchroms[allchroms != 0] # remove the dummy chromosome
for(chr in allchroms){
  PlotChromBT(pred, chr, "unadjusted")
  PlotChromBT(pred, chr, "adjusted")
}


info <- info(pred)
write.table(info, paste0(name,"/results/BubbleTree/",name,".Bubbletree.txt"), col.names=F, row.names=F, quote=F, sep="\t")

cat("\nPurity/Ploidy: ", info, "\n")

# combine multiple samples like this

# CalcBubbletreeData <- function(name){
#   loh = read.delim(paste0(name,"/results/LOH/",name,".VariantsForLOH.txt"))
#   snp.gr <- GenomicRanges::GRanges(loh$Chrom, IRanges(loh$Pos, loh$Pos), freq=loh[,"Tumor_Freq"])
#   
#   if (cnv_method=="Copywriter")
#   {
#     cnv = read.delim(paste0(name,"/results/Copywriter/",name,".seg.dat.fn"))
#   } else if (cnv_method=="HMMCopy") {
#     cnv = read.delim(paste0(name,"/results/HMMCopy/",name,".HMMCopy.20000.segments.txt.fn"))
#   }
#   
#   cnv.gr <- GenomicRanges::GRanges(cnv$Chromosome, IRanges(cnv$Start, cnv$End), num.mark=cnv$Num_Probes, seg.mean=cnv$Segment_Mean)
#   
#   r <- new("RBD", unimodal.kurtosis=-0.1)
#   
#   rbd=makeRBD(r, snp.gr,cnv.gr)
#   
#   btreepredictor <- new("BTreePredictor", rbd=rbd, max.ploidy=6, prev.grid=seq(0.2,3, by=0.01))
#   btreepredictor@config$min.segSize <- ifelse(max(btreepredictor@rbd$seg.size,na.rm=TRUE) < 0.4, 0.1, 0.4)
#   pred <- btpredict(btreepredictor)
#   return(pred)
# }
# 
# pred1 <- CalcBubbletreeData("hPDAC02_HD_LivMet-1")
# pred2 <- CalcBubbletreeData("hPDAC02_HD_PPT-1")
# 
# 
# btp <- new("BTreePlotter", max.ploidy=5, max.size=10, branch.col="gray50")
# 
# sample1 <- "hPDAC02_HD_LivMet-1"
# sample2 <- "hPDAC02_HD_PPT-1" 
# rbd1 <- pred1[[sample1]]@result$dist
# rbd2 <- pred2[[sample2]]@result$dist
# 
# srcSize <- 0.5
# trtSize <- 1
# minOver <- 1e7
# 
# arrows <- trackBTree(btp,
#                      rbd1,
#                      rbd2,
#                      min.srcSize=srcSize,
#                      min.trtSize=trtSize,
#                      min.overlap=minOver)
# 
# 
# z <- drawBTree(btp, rbd1)
# if(!is.null(arrows)) {
#   z <- z + arrows + ggplot2::labs(title=sprintf("%s -> %s", sample1, sample2))
# }
# print(z)

