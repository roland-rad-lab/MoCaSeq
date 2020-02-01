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
	for (resolution in c(10000,20000,50000,100000))
	{
	preCopywriteR(output.folder = "ref/GRCh38.p12",bin.size = resolution,ref.genome = "hg38")
	}
}

if (species=="GRCm38.p6")
{
	for (resolution in c(10000,20000,50000,100000))
	{
		preCopywriteR(output.folder = "ref/GRCm38.p6",bin.size = resolution,ref.genome = "mm10")
	}
}