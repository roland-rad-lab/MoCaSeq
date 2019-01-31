#!/usr/bin/Rscript

##########################################################################################
##
## CNV_RunCopywriter.R
##
## Run Copywriter on matched tumour-normal .bam-files using 20kb windows.
##
##########################################################################################

if(!require("CopywriteR")) install.packages("CopywriteR")
library("CopywriteR")

args <- commandArgs(TRUE)

name = args[1]
threads = args[2]

tumor_bam = paste(name,"/results/bam/",name,".Tumor.bam",sep="")
normal_bam = paste(name,"/results/bam/",name,".Normal.bam",sep="")
sample.control = data.frame(samples=c(normal_bam,tumor_bam),controls=c(normal_bam,normal_bam))

bp.param = SnowParam(workers = threads, type = "SOCK")

CopywriteR(sample.control = sample.control,
             destination.folder = file.path(paste(name,"/results/Copywriter/",sep="")),
             reference.folder = file.path("Genomes/GRCm38.p6/mm10_20kb"),
             bp.param = bp.param)

plotCNA(destination.folder = file.path(paste(name,"/results/Copywriter/",sep="")))