#!/usr/bin/Rscript

##########################################################################################
##
## LOH_MakePlots.R
##
## Plot raw data for LOH visualisation.
##
##########################################################################################

args = commandArgs(TRUE)

name=args[1] #used for naming in- and output files
repository_dir=args[2]  #location of repository

source(paste(repository_dir,"/all_GeneratePlots.R",sep=""))

setwd(paste(name,"/results/LOH",sep=""))

system(paste("mkdir ",name,"_Chromosomes",sep=""))

chrom.sizes = DefineChromSizes("Mouse")

data=paste(name,".VariantsForLOH.txt",sep="")

LOHDat = ProcessCountData(data,chrom.sizes,"LOH")
plotGlobalRatioProfile(cn=LOHDat[[1]],ChromBorders=LOHDat[[2]],cnSeg="",samplename=name,method="LOH",toolname="LOH",normalization="",y_axis="LOH",Transparency=70, Cex=0.3,outformat="pdf")

for ( i in 1:19)
{
     plotChromosomalRatioProfile(cn=LOHDat[[4]],chrom.sizes,cnSeg="",samplename=name,chromosome=i,method="LOH",toolname="LOH",SliceStart="",SliceStop="",Transparency=70, Cex=0.7, outformat="pdf")
}