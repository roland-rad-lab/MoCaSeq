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
types =  args[6]

library("CopywriteR")

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

if (species == "Human") {
	reference_files = paste(genome_dir,"/hg38_20kb",sep="")
} else if (species == "Mouse") {
	reference_files = paste(genome_dir,"/mm10_20kb",sep="")
}

bp.param = SnowParam(workers = threads, type = "SOCK")

CopywriteR(sample.control = sample.control,
             destination.folder = file.path(paste(name,"/results/Copywriter/",sep="")),
             reference.folder = file.path(reference_files),
             bp.param = bp.param)

plotCNA(destination.folder = file.path(paste(name,"/results/Copywriter/",sep="")))