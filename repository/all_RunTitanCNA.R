#!/usr/bin/Rscript

##########################################################################################
##
## all_RunTitanCNA.R
##
## Prepares all files for running Titan.
##
##########################################################################################

args = commandArgs(TRUE)

name <- args[1]
species <- args[2]
repository_dir <- args[3]
resolution <- args[4]
map_file <- args[5]
gc_file <- args[6]
exons_file <- args[7]
sequencing_type <- args[8]

suppressMessages(library(TitanCNA))
suppressMessages(library(HMMcopy))

# read in wig files and correct for GC and mappability bias
tumWig=paste(name,"/results/HMMCopy/",name,".Tumor.",resolution,".wig",sep="")
normWig=paste(name,"/results/HMMCopy/",name,".Normal.",resolution,".wig",sep="")
variants=paste(name,"/results/LOH/",name,".VariantsForLOH.txt",sep="")

if (sequencing_type == "WGS") {
	cnData=correctReadDepth(tumWig, normWig, gc_file,  map_file, genomeStyle = "NCBI") 
} else if (sequencing_type == "WES") {
	exons_file=read.delim(exons_file)
	colnames(exons_file)=c("chr","start position","stop position")
	cnData=correctReadDepth(tumWig, normWig, gc_file,  map_file, genomeStyle = "NCBI", exons_file) 
}

cnData[is.na(cnData$logR),"logR"]=0

write.table(cnData,paste(name,"/results/Titan/",name,".cnFile.txt",sep=""),sep="\t", quote=F, row.names=F,col.names=T)

variants = read.delim(variants)
variants = variants[variants[,"Normal_Freq"] <= 0.7 & variants[,"Normal_Freq"] >= 0.3,]
data = variants[,c("Chrom","Pos","Ref","Tumor_RefCount","Alt","Tumor_AltCount")]

write.table(data,paste(name,"/results/Titan/",name,".hetFile.txt",sep=""),sep="\t", quote=F, row.names=F,col.names=T)