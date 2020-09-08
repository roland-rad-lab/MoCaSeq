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

suppressPackageStartupMessages(library(BubbleTree))
suppressPackageStartupMessages(library(GenomicRanges))

loh = read.delim(paste0(name,"/results/LOH/",name,".VariantsForLOH.txt"))
snp.gr <- GenomicRanges::GRanges(loh$Chrom, IRanges(loh$Pos, loh$Pos), freq=loh[,"Tumor_Freq"])

if (cnv_method=="Copywriter")
{
	cnv = read.delim(paste0(name,"/results/Copywriter/",name,".seg.dat.fn"))
} else if (cnv_method=="HMMCopy") {
	cnv = read.delim(paste0(name,"/results/HMMCopy/",name,".HMMCopy.20000.segments.txt.fn"))
}

cnv.gr <- GenomicRanges::GRanges(cnv$Chromosome, IRanges(cnv$Start, cnv$End), num.mark=cnv$Num_Probes, seg.mean=cnv$Segment_Mean)

r <- new("RBD", unimodal.kurtosis=-0.1)

rbd=makeRBD(r, snp.gr,cnv.gr)

btreepredictor <- new("BTreePredictor", rbd=rbd, max.ploidy=6, prev.grid=seq(0.2,3, by=0.01))
btreepredictor@config$min.segSize <- ifelse(max(btreepredictor@rbd$seg.size,na.rm=TRUE) < 0.4, 0.1, 0.4)
pred <- btpredict(btreepredictor)

btreeplotter <- new("BTreePlotter", max.ploidy=10, max.size=10)
btree <- drawBTree(btreeplotter, pred@rbd)
pdf(paste0(name,"/results/BubbleTree/",name,".Bubbletree.unadjusted.pdf"))
print(btree)
garbage <- dev.off()
  
btreeplotter <- new("BTreePlotter", branch.col="gray50")
btree <- drawBTree(btreeplotter, pred@rbd.adj)
pdf(paste0(name,"/results/BubbleTree/",name,".Bubbletree.adjusted.pdf"))
print(btree)
garbage <- dev.off()

info <- info(pred)
write.table(info, paste0(name,"/results/BubbleTree/",name,".Bubbletree.txt"), col.names=F, row.names=F, quote=F, sep="\t")

cat("\nPurity/Ploidy: ", info, "\n")