#!/usr/bin/Rscript

##########################################################################################
##
## CNV_PlotCNVKit.R
##
## Plot raw data from CNVKit
##
##########################################################################################

args = commandArgs(TRUE)

name=args[1] #used for naming in- and output files
species=args[2]
repository_dir=args[3]  #location of repository
types=args[4]

source(paste(repository_dir,"/all_GeneratePlots.R",sep=""))

baseWD <- getwd()

# set the first directory
if(types == "Tumor Normal"){
  setwd(paste(name,"/results/CNVKit/matched",sep=""))
} else {
  setwd(paste(name,"/results/CNVKit/single",sep=""))
}

chrom.sizes = DefineChromSizes(species)
if (species=="Human"){
	chromosomes=22
} else if (species=="Mouse"){
	chromosomes=19
}
  
normalization=""

PlotFunction <- function(SegmentsFile, CountsFile, SampleName){
  for (y_axis in c("CNV_5","CNV_2")){
    Segments = ProcessSegmentData(segmentdata=SegmentsFile,chrom.sizes,method="CNVKit")
    Counts = ProcessCountData(countdata=CountsFile,chrom.sizes,method="CNVKit")
    
    plotGlobalRatioProfile(cn=Counts[[1]],ChromBorders=Counts[[2]],cnSeg=Segments[[1]],samplename=SampleName,method="CNV",toolname="CNVKit",normalization=normalization,y_axis=y_axis,Transparency=30, Cex=0.3,outformat="pdf")

    for (i in 1:chromosomes){
     plotChromosomalRatioProfile(cn=Counts[[4]],chrom.sizes,cnSeg=Segments[[2]],samplename=SampleName,chromosome=i,method="CNV",toolname="CNVKit",normalization=normalization,y_axis=y_axis,SliceStart="",SliceStop="",Transparency=50, Cex=0.7, outformat="pdf")
    }

    command <- paste("pdfunite ",SampleName,"_Chromosomes/",SampleName,".Chr?.CNV.CNVKit.",gsub("CNV_","",y_axis),".pdf ",SampleName,"_Chromosomes/",SampleName,".Chr??.CNV.CNVKit.",gsub("CNV_","",y_axis),".pdf ",SampleName,".Chromosomes.CNV.CNVKit.",gsub("CNV_","",y_axis),".pdf",sep="")
    system(command)
  }
}


if(types == "Tumor Normal"){
  
  system(paste("mkdir -p ",name,"_Chromosomes",sep=""))
  
  # first run the matched sample
  SegmentsFile = paste(name,".cns",sep="")
  CountsFile = paste(name,".cnr",sep="")
  PlotFunction(SegmentsFile, CountsFile, name)
      
  # second run the invididual samples again
  setwd(paste(baseWD, "/",name,"/results/CNVKit/single",sep=""))
  
  SegmentsFile = paste(name,".Tumor.cns",sep="")
  CountsFile = paste(name,".Tumor.cnr",sep="")
  subname <- paste0(name, ".Tumor")
  system(paste("mkdir -p ",subname,"_Chromosomes",sep=""))
  PlotFunction(SegmentsFile, CountsFile, subname)
  
  SegmentsFile = paste(name,".Normal.cns",sep="")
  CountsFile = paste(name,".Normal.cnr",sep="")
  subname <- paste0(name, ".Normal")
  system(paste("mkdir -p ",subname,"_Chromosomes",sep=""))
  PlotFunction(SegmentsFile, CountsFile, subname)
  
} else if(types == "Tumor"){
  system(paste("mkdir -p ",name,"_Chromosomes",sep=""))
  SegmentsFile = paste(name,".Tumor.cns",sep="")
  CountsFile = paste(name,".Tumor.cnr",sep="")
  subname <- paste0(name, ".Tumor")
  system(paste("mkdir -p ",subname,"_Chromosomes",sep=""))
  PlotFunction(SegmentsFile, CountsFile, paste0(name, ".Tumor"))
} else if(types == "Normal"){
  system(paste("mkdir -p ",name,"_Chromosomes",sep=""))
  SegmentsFile = paste(name,".Normal.cns",sep="")
  CountsFile = paste(name,".Normal.cnr",sep="")
  subname <- paste0(name, ".Normal")
  system(paste("mkdir -p ",subname,"_Chromosomes",sep=""))
  PlotFunction(SegmentsFile, CountsFile, subname)
}












