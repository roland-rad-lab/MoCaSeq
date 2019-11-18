#!/usr/bin/Rscript

##########################################################################################
##
## CNV_PlotHMMCopy.R
##
## Plot raw data from HMMCopy.
##
##########################################################################################

args = commandArgs(TRUE)

name <- args[1]
species <- args[2]
repository_dir <- args[3]
sequencing_type <- args[4]
resolution <- args[5]
map_file <- args[6]
gc_file <- args[7]
centromere_file <- args[8]
varregions_file <- args[9]

library(HMMcopy)
library(DNAcopy)
library(GenomeInfoDb)
library(naturalsort)
library(GenomicRanges)

# choose correct map_file and gc_file
map_file=gsub("20000",resolution,map_file)
gc_file=gsub("20000",resolution,gc_file)

# read in wig files and correct for GC and mappability bias
normal <- wigsToRangedData(paste(name,"/results/HMMCopy/",name,".Normal.",resolution,".wig",sep=""),gc_file,map_file)
normal$reads <- normal$reads+1
normal <- as.data.frame(correctReadcount(normal))
normal_copy=GRanges(normal$chr, IRanges(normal$start, normal$end),copy=normal$copy)

tumor <- wigsToRangedData(paste(name,"/results/HMMCopy/",name,".Tumor.",resolution,".wig",sep=""),gc_file,map_file)
tumor$reads <- tumor$reads+1
tumor <- as.data.frame(correctReadcount(tumor))
tumor_copy=GRanges(tumor$chr, IRanges(tumor$start, tumor$end),copy=tumor$copy)

# remove regions with increased variability for mice and centromere regions for humams
if (species == "Human")
{
	filtering=read.delim(centromere_file)
	flankLength=5000000
}
if (species == "Mouse")
{
	filtering=read.delim(varregions_file)
	flankLength=0
}

colnames(filtering)[1:3] <- c("space","start","end")
filtering$start <- filtering$start - flankLength
filtering$end <- filtering$end + flankLength
filtering=GRanges(filtering$space, IRanges(filtering$start, filtering$end))

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
somatic_copy$copy <- tumor_copy$copy - normal_copy$copy
somatic_tab <- as.data.frame(somatic_copy)
colnames(somatic_tab) <- c("Chrom", "Start", "End", "width", "strand", "log2Ratio")
somatic_tab <- somatic_tab[,c("Chrom", "Start", "End","log2Ratio")]
write.table(somatic_tab,paste(name,"/results/HMMCopy/",name,".HMMCopy.",resolution,".log2RR.txt", sep=""),quote=F,row.names=F,col.names=T,sep='\t')

# segmentation of the CN plot
somatic_CNA <- smooth.CNA(CNA(genomdat=somatic_tab$log2Ratio,chrom=somatic_tab$Chrom,maploc=somatic_tab$Start,data.type='logratio'))

if (sequencing_type == "lcWGS")
{
	cnv_segments <- segment(somatic_CNA,alpha=0.0001,min.width=3,undo.splits='sdundo',undo.SD=1.5,verbose=2)$output
}
if (sequencing_type == "WES" | sequencing_type == "WGS")
{
	cnv_segments <- segment(somatic_CNA,alpha=0.0001,min.width=5,undo.splits='sdundo',undo.SD=2,verbose=2)$output
}

colnames(cnv_segments) <- c("ID", "Chrom", "Start", "End", "num.mark", "Mean")
cnv_segments <- cnv_segments[,c("Chrom", "Start", "End","Mean")]
cnv_segments <- cnv_segments[naturalorder(cnv_segments$Chrom),]
write.table(cnv_segments,paste(name,"/results/HMMCopy/",name,".HMMCopy.",resolution,".segments.txt", sep=""),quote=F,row.names=F,col.names=T,sep='\t')

#start plotting
source(paste(repository_dir,"/all_GeneratePlots.R",sep=""))

setwd(paste(name,"/results/HMMCopy",sep=""))

system(paste("mkdir -p ",name,"_Chromosomes",sep=""))

chrom.sizes = DefineChromSizes(species)

if (species=="Human")
{
	chromosomes=22
} else if(species=="Mouse")
{
	chromosomes=19
}

#define normalization mode, choose from "Mode" or "MAD"

for (y_axis in c("CNV_5","CNV_2"))
{
	Segments = paste(name,".HMMCopy.",resolution,".segments.txt",sep="")
	Counts = paste(name,".HMMCopy.",resolution,".log2RR.txt",sep="")

	Segments = ProcessSegmentData(segmentdata=Segments,chrom.sizes,method="HMMCopy")
	Counts = ProcessCountData(countdata=Counts,chrom.sizes,method="HMMCopy")

	plotGlobalRatioProfile(cn=Counts[[1]],ChromBorders=Counts[[2]],cnSeg=Segments[[1]],samplename=name,method="CNV",toolname="HMMCopy",normalization="",y_axis=y_axis,Transparency=30, Cex=0.3,outformat="pdf")

	for (i in 1:chromosomes)
	{
    plotChromosomalRatioProfile(cn=Counts[[4]],chrom.sizes,cnSeg=Segments[[2]],samplename=name,chromosome=i,method="CNV",toolname="HMMCopy",normalization="",y_axis=y_axis,SliceStart="",SliceStop="",Transparency=50, Cex=0.7, outformat="pdf")
	}
	
	system(paste("pdfunite ",name,"_Chromosomes/",name,".Chr?.CNV.HMMCopy.",gsub("CNV_","",y_axis),".pdf ",name,"_Chromosomes/",name,".Chr??.CNV.HMMCopy.",gsub("CNV_","",y_axis),".pdf ",name,".Chromosomes.CNV.HMMCopy.",gsub("CNV_","",y_axis),".pdf",sep=""))
}