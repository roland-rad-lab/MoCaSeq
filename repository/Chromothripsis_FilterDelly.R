#!/usr/bin/Rscript

##########################################################################################
##
## Chromothripsis_FilterDelly.R
##
## Takes an annotated DELLY output file and performs filtering steps.
##
##########################################################################################

options(warn=-1)
if(!require("optparse")) install.packages("optparse")
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-n", "--name"),type="character",default=NULL,help="name"),
  make_option(c("-i", "--input"),type="character",default=NULL,help="annotated rearrangement list"),
  make_option(c("-o", "--intraChromOnly"),type="logical",default=T,help="logical flag; If set to TRUE, only intrachromosomal rearrangements are kept."),
  make_option(c("-c", "--coverageFilter"),type="numeric",default=NULL,help="estimated coverage filter")
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
FilterInterChromosomalRearrangements = opt$intraChromOnly
CoverageFilter = opt$coverageFilter
name = opt$name

# main program
print(paste0("Coverage Filter set to ",CoverageFilter))
if(FilterInterChromosomalRearrangements) print(paste("Filter only intrachromosomal Rearrangements"))
print("Read Input File")

# apply filtering steps
if(FilterInterChromosomalRearrangements) input_list = input_list[input_list$donChr==input_list$accChr,]
input_list <- input_list[input_list$control.varFreq == 0 & input_list$tumor.varFreq >= 0.2 & input_list$distances >= 6000,]

### Distance set because of huge FragmentSizes
### Control coverage set as filter for repeat regions , changed to 40 x
ImpList = input_list[input_list$Precision == "IMPRECISE" & input_list$distances > 1000000 & input_list$control.coverage <= CoverageFilter &
          input_list$MAPQ >= 30,]
PreList = input_list[input_list$Precision == "PRECISE",]
input_list = rbind(ImpList,PreList)
input_list = unique(input_list[order(input_list[,1],input_list[,2]),])

print("Writing Output")
write.table(input_list,paste(name,"/results/Delly/",name,".breakpoints.filtered.tab",sep=""),quote=F,row.names=F,col.names=T,sep='\t')