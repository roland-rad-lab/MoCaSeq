#!/usr/bin/Rscript

##########################################################################################
##
## CNV_CopywriterGetModeCorrectionFactor.R
##
## Shift datapoints and segments from CopywriteR using the calculated correction factor.
##
##########################################################################################

args <- commandArgs(TRUE)

name = args[1]
runmode = args[2]
type = args[3]

path = getwd()

### Descend into Copywriter results folder
setwd(paste(path,"/",name,"/results/Copywriter/",sep=""))
ShiftData = read.table(paste(name,".KDEestimate.txt",sep=""),header=T,sep="\t")
segmentData = read.table(paste(name,".Copywriter.segments.MAD.txt",sep=""),header=T,sep="\t")

Shift = ShiftData$Shift
segmentData$seg.mean = segmentData$Mean-Shift
segmentData=segmentData[,c("Chrom", "Start", "End", "Mean")]
write.table(segmentData,file=paste(name,".Copywriter.segments.Mode.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)

logReadCounts = read.table("CNAprofiles/log2_read_counts.igv",header=T,sep="\t")

if (runmode == "MS") {
	TumorLogReadCounts = grep("Tumor",colnames(logReadCounts),value=T)
	NormalLogReadCounts = grep("Normal",colnames(logReadCounts),value=T)
	logReadCounts$Copy = logReadCounts[,TumorLogReadCounts]-logReadCounts[,NormalLogReadCounts]
	logReadCountsMode = logReadCounts

	logReadCountsMode[,TumorLogReadCounts] = logReadCountsMode[,TumorLogReadCounts]-Shift
	logReadCountsMode$Copy = logReadCountsMode[,TumorLogReadCounts]-logReadCountsMode[,NormalLogReadCounts]
	logReadCountsMode = logReadCountsMode[,c("Chromosome","Start","End","Copy")]
	colnames(logReadCountsMode) = c("Chrom","Start","End","log2Ratio")
	write.table(logReadCountsMode,file=paste(name,".Copywriter.log2RR.Mode.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)

	logReadCounts = logReadCounts[,c("Chromosome","Start","End","Copy")]
	colnames(logReadCounts) = c("Chrom","Start","End","log2Ratio")
	write.table(logReadCounts,file=paste(name,".Copywriter.log2RR.MAD.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
} else if (runmode == "SS") {
	Selection = paste("log2.",gsub("-",".",name),".",type,".bam",sep="")
	
	logReadCounts$Copy = logReadCounts[,Selection]
	logReadCountsMode = logReadCounts
	logReadCountsMode[,Selection] = logReadCountsMode[,Selection]-Shift
	logReadCountsMode = logReadCountsMode[,c("Chromosome","Start","End","Copy")]
	colnames(logReadCountsMode) = c("Chrom","Start","End","log2Ratio")
	write.table(logReadCountsMode,file=paste(name,".Copywriter.log2RR.Mode.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)

	logReadCounts = logReadCounts[,c("Chromosome","Start","End","Copy")]
	colnames(logReadCounts) = c("Chrom","Start","End","log2Ratio")
	write.table(logReadCounts,file=paste(name,".Copywriter.log2RR.MAD.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
}

setwd(path)