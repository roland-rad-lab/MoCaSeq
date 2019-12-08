#!/usr/bin/Rscript

##########################################################################################
##
## CNV_PlotaCGH.R
##
## Plot raw data from Agilent Genomic Workbench.
##
##########################################################################################

args = commandArgs(TRUE)

name=args[1] #used for naming in- and output files
species=args[2]
repository_dir=args[3]  #location of repository

source(paste(repository_dir,"/all_GeneratePlots.R",sep=""))

system(paste("mkdir -p ",name,"_Chromosomes",sep=""))

chrom.sizes = DefineChromSizes(species)

if (species=="Human")
{
	chromosomes=22
} else if (species=="Mouse")
{
	chromosomes=19
}

#define normalization mode, choose from "Mode" or "MAD"
normalization="Mode"

for (y_axis in c("CNV_5","CNV_2"))
{
	Segments = paste(name,"-PPT_244k_flat.txt",sep="")	
	Counts = paste(name,"-PPTforValidation.txt",sep="")

	Counts = ProcessCountData(countdata=Counts,chrom.sizes,method="aCGH")
	ProcessedSegments = PreProcess_aCGH(Segments,chrom.sizes,Counts[[3]])
	Segments = ProcessSegmentData(segmentdata=ProcessedSegments,chrom.sizes,method="aCGH")
	
	plotGlobalRatioProfile(cn=Counts[[1]],ChromBorders=Counts[[2]],cnSeg=Segments[[1]],samplename=name,method="CNV",toolname="aCGH",normalization=normalization,y_axis=y_axis,Transparency=30, Cex=0.3,outformat="pdf")

	for (i in 1:chromosomes)
	{
    plotChromosomalRatioProfile(cn=Counts[[4]],chrom.sizes,cnSeg=Segments[[2]],samplename=name,chromosome=i,method="CNV",toolname="aCGH",normalization=normalization,y_axis=y_axis,SliceStart="",SliceStop="",Transparency=50, Cex=0.7, outformat="pdf")
	}

	system(paste("pdfunite ",name,"_Chromosomes/",name,".Chr?.CNV.aCGH.Mode.",gsub("CNV_","",y_axis),".pdf ",name,"_Chromosomes/",name,".Chr??.CNV.aCGH.Mode.",gsub("CNV_","",y_axis),".pdf ",name,".Chromosomes.CNV.aCGH.Mode.",gsub("CNV_","",y_axis),".pdf",sep=""))
}