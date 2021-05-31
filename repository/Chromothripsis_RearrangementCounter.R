#!/usr/bin/Rscript

##########################################################################################
##
## Chromothripsis_RearrangementCounter.R
##
## Counts the number of rearrangements for the provided chromosome.
##
##########################################################################################

options(warn=-1)
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-i", "--input"),type="character",default=NULL,help="rearrangement list"),
  make_option(c("-c", "--chrom"),type="character",default=NULL,help="chromosome to analyze")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input_list <- read.csv(opt$input,header=T,sep="\t",stringsAsFactors=F)
input_list = input_list[input_list$donChr==opt$chrom & input_list$accChr==opt$chrom,]
cat(nrow(input_list))