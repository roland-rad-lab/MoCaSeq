#!/usr/bin/Rscript

##########################################################################################
##
## CNV_CopywriterGetRawData.R
##
## Extract datapoints and segments from the Rdata object provided by CopywriteR.
##
##########################################################################################

args <- commandArgs(TRUE)

name = args[1]
runmode = args[2]
type = args[3]

# Keep track of working directory
path = getwd()

# Descend into Copywriter results folder
setwd(paste(path,"/",name,"/results/Copywriter/CNAprofiles/",sep=""))

#load raw segment data from file, change around coluns and export them
load("segment.Rdata")
segmentData = segment.CNA.object$output

if (runmode == "MS") {
	Selection = unique(grep("Normal",grep("Tumor",segmentData$ID,value=T),value=T))
} else if (runmode == "SS") {
	Selection = paste("log2.",gsub("-",".",name),".",type,".bam.vs.none",sep="")
} 

segmentData = segmentData[segmentData$ID==Selection,]
chromosomes = unique(segmentData$chrom)
# Exclude sex chromosomes
chromosomes = chromosomes[chromosomes<=22]
for(chromosome in chromosomes)
   {
   segmentMean = as.numeric(segmentData[segmentData$chrom==chromosome,"seg.mean"])
   write.table(segmentMean,file=paste("../",name,".SegmentsChromosome",chromosome,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
   }
segmentData = segmentData[,c("chrom","loc.start","loc.end","seg.mean")]
colnames(segmentData)=c("Chrom","Start","End","Mean")
segmentData$Start = floor(segmentData$Start)
write.table(segmentData,file=paste("../",name,".Copywriter.segments.MAD.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
setwd(path)