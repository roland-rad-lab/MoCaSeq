#!/usr/bin/Rscript

##########################################################################################
## 
## Chromothripsis_AnnotateRatios.R
##
## Explicity annotate read support (pair/split read) for all rearrangements detected by Delly.
##
##########################################################################################

options(warn=-1)
if(!require("optparse")) install.packages("optparse")
library(optparse)

option_list = list(
  make_option(c("-i", "--input"),type="character",default=NULL,help="delly rearrangement vcf file")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if(is.null(opt$input)){
  print_help(opt_parser)
  stop("You have to specify a rearrangement list.", call.=FALSE)
} else{
  tryCatch({
    input_list <- read.csv(opt$input,header=T,sep="\t",stringsAsFactors=F)
  },
  error = function(e){
    stop("Error reading rearrangement list.")
  })
}


# main program
control.pairs.ratio <- input_list$control.DV / (input_list$control.DV + input_list$control.DR)
control.pairs.ratio[is.na(control.pairs.ratio)] = 0
control.SR.ratio <- input_list$control.RV / (input_list$control.RV + input_list$control.RR)
control.SR.ratio[is.na(control.SR.ratio)] = 0
tumor.pairs.ratio <- input_list$tumor.DV / (input_list$tumor.DV + input_list$tumor.DR)
tumor.pairs.ratio[is.na(tumor.pairs.ratio)] = 0
tumor.SR.ratio <- input_list$tumor.RV / (input_list$tumor.RV + input_list$tumor.RR)
tumor.SR.ratio[is.na(tumor.SR.ratio)] = 0
control.coverage <- input_list$control.DR + input_list$control.DV + input_list$control.RR + input_list$control.RV
tumor.coverage <- input_list$tumor.DR + input_list$tumor.DV + input_list$tumor.RR + input_list$tumor.RV
control.varFreq <- (input_list$control.DV + input_list$control.RV) / control.coverage
control.varFreq[is.na(control.varFreq)] = 0
tumor.varFreq <- (input_list$tumor.DV + input_list$tumor.RV) / tumor.coverage
tumor.varFreq[is.na(tumor.varFreq)] = 0
distances <- input_list$accPos - input_list$donPos

output_list <- cbind(input_list,control.pairs.ratio,control.SR.ratio,tumor.pairs.ratio,tumor.SR.ratio,control.varFreq,tumor.varFreq,control.coverage,tumor.coverage,distances)

write.table(output_list,"",quote=F,row.names=F,col.names=T,sep='\t')