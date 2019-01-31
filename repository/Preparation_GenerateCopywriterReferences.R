#!/usr/bin/Rscript

##########################################################################################
##
## Preparation_GenerateCopywriterReferences.R
##
## Creates Copywriter reference files for 20kb windows.
##
##########################################################################################

if(!require("CopywriteR")) install.packages("CopywriteR")
library("CopywriteR")

preCopywriteR(output.folder= "Genomes/GRCm38.p6",bin.size = 20000,ref.genome = "mm10")
