#!/usr/bin/Rscript

##########################################################################################
##
## CNV_RunCopywriter.R
##
## Run Copywriter on matched tumour-normal .bam-files using 20kb windows.
##
##########################################################################################

args <- commandArgs(TRUE)

name = args[1]
species = args[2]
threads = args[3]
runmode = args[4]
genome_dir = args[5]
centromere_file <- args[6]
varregions_file <- args[7]
resolution <- args[8]
types <- args[9]

if (resolution=="NULL") { resolution=20000 } 

suppressPackageStartupMessages(library(CopywriteR))
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(naturalsort))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))

tumor_bam = paste(name,"/results/bam/",name,".Tumor.bam",sep="")
normal_bam = paste(name,"/results/bam/",name,".Normal.bam",sep="")

if (runmode == "MS") {
	sample.control = data.frame(samples=c(normal_bam,tumor_bam),controls=c(normal_bam,normal_bam))
}

if (runmode == "SS") {
	if (types == "Tumor") {
		sample.control = data.frame(samples=c(tumor_bam),controls=c(tumor_bam))
	} else if (types == "Normal") {
		sample.control = data.frame(samples=c(normal_bam),controls=c(normal_bam))
	}
}

resolution=as.numeric(as.character(resolution))/1000
destination.folder = file.path(paste(name,"/results/Copywriter/",sep=""))

# Copywriter will break if the destination folder already has the folders from a previous run (this is only an issue for rerunning Copywriter)
if(dir.exists(paste0(destination.folder, "/CNAprofiles"))){
  unlink(paste0(destination.folder, "/CNAprofiles/"), recursive = T)
}

if (species == "Human") {
	reference_files = paste(genome_dir,"/hg38_",resolution,"kb",sep="")
} else if (species == "Mouse") {
	reference_files = paste(genome_dir,"/mm10_",resolution,"kb",sep="")
}

bp.param = SnowParam(workers = threads, type = "SOCK")

CopywriteR(sample.control = sample.control,
             destination.folder = destination.folder,
             reference.folder = file.path(reference_files),
             bp.param = bp.param)

# filter counts
log2.reads=read.table(paste0(name,"/results/Copywriter/CNAprofiles/log2_read_counts.igv"), header=T, sep="\t",check.names=FALSE)

if (runmode == "SS") {
  if (types == "Tumor") {
    log2.reads.GR=GRanges(log2.reads$Chromosome, IRanges(log2.reads$Start, log2.reads$End),Feature=as.character(log2.reads$Feature), Normal=NA,Tumor=log2.reads[,5])
  } else if (types == "Normal") {
    log2.reads.GR=GRanges(log2.reads$Chromosome, IRanges(log2.reads$Start, log2.reads$End),Feature=as.character(log2.reads$Feature), Normal=log2.reads[,5],Tumor=NA)
  }
} else {
  log2.reads.GR=GRanges(log2.reads$Chromosome, IRanges(log2.reads$Start, log2.reads$End),Feature=as.character(log2.reads$Feature), Normal=log2.reads[,5],Tumor=log2.reads[,6])
}

file.copy(paste0(name,"/results/Copywriter/CNAprofiles/log2_read_counts.igv"),paste0(name,"/results/Copywriter/CNAprofiles/log2_read_counts_raw.igv"), overwrite=T)  # backup copy

# remove regions with increased variability for mice and centromere regions for humams
if (species == "Human"){
	filter=read.delim(centromere_file)
	flankLength=5000000
}
if (species == "Mouse"){
	filter=read.delim(varregions_file)
	flankLength=0
}

colnames(filter)[1:3] <- c("Chromosome","Start","End")
filter$Start <- filter$Start - flankLength
filter$End <- filter$End + flankLength
filter=GRanges(filter$Chromosome, IRanges(filter$Start, filter$End))

hits <- findOverlaps(query = log2.reads.GR, subject = filter)
ind <- queryHits(hits)
message("Removed ", length(ind), " bins near centromeres (human) or variable regions (mouse).")

# remove those regions (ind = 0 would remove all)
if(length(ind) != 0){
  log2.reads.GR=(log2.reads.GR[-ind, ])
}

log2.reads.fixed=as.data.frame(log2.reads.GR)
log2.reads.fixed=log2.reads.fixed[,c("seqnames", "start", "end", "Feature", "Normal","Tumor")]
colnames(log2.reads.fixed)=c("Chromosome", "Start", "End", "Feature",paste0("log2.",name,".Normal.bam"),paste0("log2.",name,".Tumor.bam"))

write.table(log2.reads.fixed,paste0(name,"/results/Copywriter/CNAprofiles/log2_read_counts.igv"), sep="\t", quote=F, row.names=F,col.names=T)

# remove the hardcoded path (/var/pipeline/) in the input data (else noone can repeat this analysis outside of docker)
outfolder <- paste0(destination.folder, "/CNAprofiles")

load(file.path(outfolder, "input.Rdata"))
samplesDT <- data.table(inputStructure$sample.control)
samplesDT[, samples := gsub("/var/pipeline/", "", samples)]
samplesDT[, controls := gsub("/var/pipeline/", "", controls)]
inputStructure$sample.control <- data.frame(samplesDT)
save(inputStructure, file=paste0(outfolder, "/input.Rdata"))

# Plot
plotCNA(destination.folder = file.path(paste(name,"/results/Copywriter/",sep=""))) # will also create segRdataFile.Rdata

# rename chromosomes in segment.Rdata
segRdataFile <- paste0(name,"/results/Copywriter/CNAprofiles/segment.Rdata")
load(segRdataFile)
file.copy(segRdataFile,paste0(name,"/results/Copywriter/CNAprofiles/segment_raw.Rdata"), overwrite=T) # backup copy
segmentData = segment.CNA.object$output

if (species == "Human") {
  segmentData$chrom[segmentData$chrom == 23] <- "X"
  segmentData$chrom[segmentData$chrom == 24] <- "Y"
} else if (species == "Mouse") {
  segmentData$chrom[segmentData$chrom == 20] <- "X"
  segmentData$chrom[segmentData$chrom == 21] <- "Y"
}

segment.CNA.object$output <- segmentData
save(segment.CNA.object, file=segRdataFile)

