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
  
suppressPackageStartupMessages(library(TitanCNA))
suppressPackageStartupMessages(library(HMMcopy))
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(GenomicRanges))

correctReadDepth = function (tumWig, normWig, gcWig, mapWig, genomeStyle = "NCBI", 
    targetedSequence = NULL) 
{
    message("Reading GC and mappability files")
    gc <- wigToGRanges(gcWig)
    map <- wigToGRanges(mapWig)
    message("Loading tumour file:", tumWig)
    tumour_reads <- wigToGRanges(tumWig)
    message("Loading normal file:", normWig)
    normal_reads <- wigToGRanges(normWig)
    seqlevelsStyle(gc) <- genomeStyle
    seqlevelsStyle(map) <- genomeStyle
    seqlevelsStyle(tumour_reads) <- genomeStyle
    seqlevelsStyle(normal_reads) <- genomeStyle
    gc <- gc[seqnames(gc) %in% seqnames(tumour_reads)]
    map <- map[seqnames(map) %in% seqnames(tumour_reads)]
    samplesize <- 50000
    if (!is.null(targetedSequence)) {
        message("Analyzing targeted regions...")
        targetIR <- GRanges(ranges = IRanges(start = targetedSequence[, 
            2], end = targetedSequence[, 3]), seqnames = targetedSequence[, 
            1])
        names(targetIR) <- setGenomeStyle(seqlevels(targetIR), genomeStyle)
        hits <- findOverlaps(query = tumour_reads, subject = targetIR)
        keepInd <- unique(queryHits(hits))
        tumour_reads <- tumour_reads[keepInd, ]
        normal_reads <- normal_reads[keepInd, ]
        gc <- gc[keepInd, ]
        map <- map[keepInd, ]
        samplesize <- min(ceiling(nrow(tumour_reads) * 0.1), 
            samplesize)
    }
    tumour_reads$gc <- gc$value
    tumour_reads$map <- map$value
    colnames(values(tumour_reads)) <- c("reads", "gc", "map")
    normal_reads$gc <- gc$value
    normal_reads$map <- map$value
    colnames(values(normal_reads)) <- c("reads", "gc", "map")
    message("Correcting Tumour")
    tumour_copy <- correctReadcount(tumour_reads, samplesize = samplesize)
    message("Correcting Normal")
    normal_copy <- correctReadcount(normal_reads, samplesize = samplesize)
    message("Normalizing Tumour by Normal")
    tumour_copy$copy <- tumour_copy$copy - normal_copy$copy
    rm(normal_copy)
    temp <- cbind(chr = as.character(seqnames(tumour_copy)), 
        start = start(tumour_copy), end = end(tumour_copy), logR = tumour_copy$copy)
    temp <- as.data.frame(temp, stringsAsFactors = FALSE)
    mode(temp$start) <- "numeric"
    mode(temp$end) <- "numeric"
    mode(temp$logR) <- "numeric"
    return(temp)
}

# read in wig files and correct for GC and mappability bias
tumWig=paste(name,"/results/HMMCopy/",name,".Tumor.",resolution,".wig",sep="")
normWig=paste(name,"/results/HMMCopy/",name,".Normal.",resolution,".wig",sep="")
variants=paste(name,"/results/LOH/",name,".VariantsForLOH.txt",sep="")

if (sequencing_type == "WGS") {
	cnData=correctReadDepth(tumWig, normWig, gc_file,  map_file, genomeStyle = "NCBI") 
} else if (sequencing_type == "WES") {
	exons_file=read.delim(exons_file)
	exons_file=exons_file[,1:3]
	colnames(exons_file)=c("chr","start","end")
	cnData=correctReadDepth(tumWig, normWig, gc_file,  map_file, genomeStyle = "NCBI", exons_file) 
}

cnData[is.na(cnData$logR),"logR"]=0

write.table(cnData,paste(name,"/results/Titan/",name,".cnFile.txt",sep=""),sep="\t", quote=F, row.names=F,col.names=T)

variants = read.delim(variants)
variants = variants[variants[,"Normal_Freq"] <= 0.7 & variants[,"Normal_Freq"] >= 0.3,]
data = variants[,c("Chrom","Pos","Ref","Tumor_RefCount","Alt","Tumor_AltCount")]

write.table(data,paste(name,"/results/Titan/",name,".hetFile.txt",sep=""),sep="\t", quote=F, row.names=F,col.names=T)