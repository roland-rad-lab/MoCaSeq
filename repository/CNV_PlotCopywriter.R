#!/usr/bin/Rscript

##########################################################################################
##
## CNV_PlotCopywriter.R
##
## Plot raw data from Copywriter.
##
##########################################################################################

args = commandArgs(TRUE)

name=args[1] #used for naming in- and output files
species=args[2]
repository_dir=args[3]  #location of repository

source(paste(repository_dir,"/all_GeneratePlots.R",sep=""))

setwd(paste(name,"/results/Copywriter",sep=""))

system(paste("mkdir -p ",name,"_Chromosomes",sep=""))

chrom.sizes = DefineChromSizes(species)

if (species=="Human"){
	chromosomes=22
} else if (species=="Mouse"){
	chromosomes=19
}

#define normalization mode, choose from "Mode" or "MAD"
normalization="Mode"

for (y_axis in c("CNV_5","CNV_2"))
{
	Segments = paste(name,".Copywriter.segments.",normalization,".txt",sep="")
	Counts = paste(name,".Copywriter.log2RR.",normalization,".txt",sep="")

	Segments = ProcessSegmentData(segmentdata=Segments,chrom.sizes,method="Copywriter")
	Counts = ProcessCountData(countdata=Counts,chrom.sizes,method="Copywriter")

	plotGlobalRatioProfile(cn=Counts[[1]],ChromBorders=Counts[[2]],cnSeg=Segments[[1]],samplename=name,method="CNV",toolname="Copywriter",normalization=normalization,y_axis=y_axis,Transparency=30, Cex=0.3,outformat="pdf")

	for (i in 1:chromosomes)
	{
    plotChromosomalRatioProfile(cn=Counts[[4]],chrom.sizes,cnSeg=Segments[[2]],samplename=name,chromosome=i,method="CNV",toolname="Copywriter",normalization=normalization,y_axis=y_axis,SliceStart="",SliceStop="",Transparency=50, Cex=0.7, outformat="pdf")
	}

	system(paste("pdfunite ",name,"_Chromosomes/",name,".Chr?.CNV.Copywriter.Mode.",gsub("CNV_","",y_axis),".pdf ",name,"_Chromosomes/",name,".Chr??.CNV.Copywriter.Mode.",gsub("CNV_","",y_axis),".pdf ",name,".Chromosomes.CNV.Copywriter.Mode.",gsub("CNV_","",y_axis),".pdf",sep=""))
}