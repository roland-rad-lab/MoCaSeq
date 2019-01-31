#!/usr/bin/Rscript

##########################################################################################
##
## CNV_PlotHMMCopy.R
##
## Plot raw data from HMMCopy.
##
##########################################################################################

if(!require("HMMcopy")) install.packages("HMMcopy")
if(!require("DNAcopy")) install.packages("DNAcopy")
if(!require("GenomeInfoDb")) install.packages("GenomeInfoDb")
library(HMMcopy)
library(DNAcopy)
library(GenomeInfoDb)

args = commandArgs(TRUE)

name <- args[1]
repository_dir <- args[2]
resolution <- args[3]

map_file=paste("Genomes/GRCm38.p6/GRCm38.p6.map.",resolution,".wig", sep="")
gc_file=paste("Genomes/GRCm38.p6/GRCm38.p6.gc.",resolution,".wig", sep="")

# read in wig files and correct for GC and mappability bias
normal_copy <- correctReadcount(wigsToRangedData(paste(name,"/results/HMMCopy/",name,".Normal.",resolution,".wig",sep=""),gc_file,map_file))
tumor_copy <- correctReadcount(wigsToRangedData(paste(name,"/results/HMMCopy/",name,".Tumor.",resolution,".wig",sep=""),gc_file,map_file))

# computation of the copy number states from the log fold change
somatic_copy <- tumor_copy
somatic_copy$copy <- tumor_copy$copy - normal_copy$copy
somatic_tab <- as.data.frame(somatic_copy)
colnames(somatic_tab) <- c("Chrom", "Start", "End", "width", "reads", "gc", "map", "valid", "ideal", "cor.gc", "cor.map", "log2Ratio")
somatic_tab <- somatic_tab[,c("Chrom", "Start", "End","log2Ratio")]
write.table(somatic_tab,paste(name,"/results/HMMCopy/",name,".HMMCopy.",resolution,".log2RR.txt", sep=""),quote=F,row.names=F,col.names=T,sep='\t')

# segmentation of the CN plot
somatic_CNA <- smooth.CNA(CNA(genomdat=somatic_tab$log2Ratio,chrom=somatic_tab$Chrom,maploc=somatic_tab$Start,data.type='logratio'))
cnv_segments <- segment(somatic_CNA,alpha=0.0001,min.width=5,undo.splits='sdundo',undo.SD=2,verbose=2)$output
colnames(cnv_segments) <- c("ID", "Chrom", "Start", "End", "num.mark", "Mean")
cnv_segments <- cnv_segments[,c("Chrom", "Start", "End","Mean")]
write.table(cnv_segments,paste(name,"/results/HMMCopy/",name,".HMMCopy.",resolution,".segments.txt", sep=""),quote=F,row.names=F,col.names=T,sep='\t')

#start plotting
source(paste(repository_dir,"/all_GeneratePlots.R",sep=""))

setwd(paste(name,"/results/HMMCopy",sep=""))

system(paste("mkdir ",name,"_Chromosomes",sep=""))

chrom.sizes = DefineChromSizes("Mouse")

#define normalization mode, choose from "Mode" or "MAD"

for (y_axis in c("CNV_5","CNV_2"))
{
	Segments = paste(name,".HMMCopy.",resolution,".segments.txt",sep="")
	Counts = paste(name,".HMMCopy.",resolution,".log2RR.txt",sep="")

	Segments = ProcessSegmentData(segmentdata=Segments,chrom.sizes,method="HMMCopy")
	Counts = ProcessCountData(countdata=Counts,chrom.sizes,method="HMMCopy")

	plotGlobalRatioProfile(cn=Counts[[1]],ChromBorders=Counts[[2]],cnSeg=Segments[[1]],samplename=name,method="CNV",toolname="HMMCopy",normalization="",y_axis=y_axis,Transparency=70, Cex=0.3,outformat="pdf")

	for ( i in 1:19)
	{
    plotChromosomalRatioProfile(cn=Counts[[4]],chrom.sizes,cnSeg=Segments[[2]],samplename=name,chromosome=i,method="CNV",toolname="HMMCopy",normalization="",y_axis=y_axis,SliceStart="",SliceStop="",Transparency=70, Cex=0.7, outformat="pdf")
	}
}