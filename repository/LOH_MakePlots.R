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
species=args[2]
repository_dir=args[3]  #location of repository

source(paste(repository_dir,"/all_GeneratePlots.R",sep=""))

setwd(paste(name,"/results/LOH",sep=""))

system(paste("mkdir -p ",name,"_Chromosomes",sep=""))

chrom.sizes = DefineChromSizes(species)

if (species=="Human")
{
	chromosomes=21
} else if (species=="Mouse")
{
	chromosomes=19
}

data=paste(name,".VariantsForLOH.txt",sep="")

LOHDat = ProcessCountData(data,chrom.sizes,"LOH")
plotGlobalRatioProfile(cn=LOHDat[[1]],ChromBorders=LOHDat[[2]],cnSeg="",samplename=name,method="LOH",toolname="LOH",normalization="",y_axis="LOH",Transparency=70, Cex=0.3,outformat="pdf")

for (i in 1:chromosomes)
{
     plotChromosomalRatioProfile(cn=LOHDat[[4]],chrom.sizes,cnSeg="",samplename=name,chromosome=i,method="LOH",toolname="LOH",SliceStart="",SliceStop="",Transparency=70, Cex=0.7, outformat="pdf")
}

system(paste("pdfunite ",name,"_Chromosomes/",name,".Chr?.LOH.LOH.pdf ",name,"_Chromosomes/",name,".Chr??.LOH.LOH.pdf ",name,".Chromosomes.LOH.LOH.pdf",sep=""))

data=paste(name,".VariantsForLOHGermline.txt",sep="")

LOH_GermlineDat = ProcessCountData(data,chrom.sizes,"LOH_Germline")
plotGlobalRatioProfile(cn=LOH_GermlineDat[[1]],ChromBorders=LOH_GermlineDat[[2]],cnSeg="",samplename=name,method="LOH_Germline",toolname="LOH_Germline",normalization="",y_axis="LOH_Germline",Transparency=70, Cex=0.3,outformat="pdf")

for (i in 1:chromosomes)
{
     plotChromosomalRatioProfile(cn=LOH_GermlineDat[[4]],chrom.sizes,cnSeg="",samplename=name,chromosome=i,method="LOH_Germline",toolname="LOH_Germline",SliceStart="",SliceStop="",Transparency=70, Cex=0.7, outformat="pdf")
}

system(paste("pdfunite ",name,"_Chromosomes/",name,".Chr?.LOH_Germline.LOH_Germline.pdf ",name,"_Chromosomes/",name,".Chr??.LOH_Germline.LOH_Germline.pdf ",name,".Chromosomes.LOH_Germline.LOH_Germline.pdf",sep=""))
