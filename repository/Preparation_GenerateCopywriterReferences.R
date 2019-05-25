#!/usr/bin/Rscript

##########################################################################################
##
## Preparation_GenerateCopywriterReferences.R
##
## Creates Copywriter reference files for 20kb windows.
##
##########################################################################################

args <- commandArgs(TRUE)

species = args[1]

library("CopywriteR")

if (species=="GRCh38.p12")
{
	preCopywriteR(output.folder = "ref/GRCh38.p12",bin.size = 20000,ref.genome = "hg38")
}

if (species=="GRCm38.p6")
{
	preCopywriteR(output.folder = "ref/GRCm38.p6",bin.size = 20000,ref.genome = "mm10")
}