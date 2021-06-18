#!/usr/bin/Rscript

##########################################################################################
## 
## Chromothripsis_AnnotateRatios.R
##
## Explicity annotate read support (pair/split read) for all rearrangements detected by Delly and perform filtering steps
##
##########################################################################################

options(warn=-1)
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

option_list = list(
  make_option(c("-n", "--name"),type="character",default=NULL,help="sample name"),
  make_option(c("-i", "--input"),type="character",default=NULL,help="rearrangement tab separated file (from Chromothripsis_FormatTable.sh)"),
  make_option(c("-o", "--intraChromOnly"),type="logical",default=T,help="logical flag; If set to TRUE, only intrachromosomal rearrangements are kept."),
  make_option(c("-c", "--coverageFilter"),type="numeric",default=NULL,help="estimated coverage filter"),
  make_option(c("-v", "--filterTumorVarFreq"),type="numeric",default=0.2,help="Variants will be filtered based on tumor variant freq >= x (default is 0.2)"),
  make_option(c("-d", "--distance"),type="numeric",default=6000,help="Variants will be filtered based on distance between breakpoints >= x (default is 6000)")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)
if (is.null(opt$name)) stop("Need --name")
if (is.null(opt$input)) stop("Need --input")
if (is.null(opt$intraChromOnly)) stop("Need --outfile.")
if (is.null(opt$coverageFilter)) stop("Need --coverageFilter")

name <- opt$name
input <- opt$input
FilterInterChromosomalRearrangements <- opt$intraChromOnly #T
CoverageFilter <- opt$coverageFilter # calculated externally
TumorVarFreq <- opt$filterTumorVarFreq # 0.2
maxDist <- opt$distance #6000

# FORMATTING
inputDT <- fread(input,header=T,sep="\t",stringsAsFactors=F)

# rename (extensively) like this, to facilitate future coding
setnames(inputDT, "GEN[Normal].DR", "control.DR")
setnames(inputDT, "GEN[Normal].CN", "control.DV")
setnames(inputDT, "GEN[Normal].RR", "control.RR")
setnames(inputDT, "GEN[Normal].RV", "control.RV")
setnames(inputDT, "GEN[Tumor].DR", "tumor.DR")
setnames(inputDT, "GEN[Tumor].CN", "tumor.DV")
setnames(inputDT, "GEN[Tumor].RR", "tumor.RR")
setnames(inputDT, "GEN[Tumor].RV", "tumor.RV")
setnames(inputDT, c("CHROM", "POS", "CHR2", "POS2"), c("donChr", "donPos", "accChr", "accPos"))
setnames(inputDT, "CT", "Type")

# set types to avoid warnings
inputDT[, accPos := as.integer(accPos)]
inputDT[, accChr := as.integer(accChr)]

# empty acceptor CHROM+POS is equal to changes on the same chromosome
inputDT[accChr == "" | is.na(accChr), accPos := END]
inputDT[accChr == "" | is.na(accChr), accChr := donChr]
inputDT[, END := NULL]

# ANNOTATION
inputDT[, control.pairs.ratio := control.DV / (control.DV + control.DR)]
inputDT[is.na(control.pairs.ratio), control.pairs.ratio := 0]
inputDT[, control.SR.ratio := control.RV / (control.RV + control.RR)]
inputDT[is.na(control.SR.ratio), control.SR.ratio := 0]

inputDT[, tumor.pairs.ratio := tumor.DV / (tumor.DV + tumor.DR)]
inputDT[is.na(tumor.pairs.ratio), tumor.pairs.ratio := 0]
inputDT[, tumor.SR.ratio := tumor.RV / (tumor.RV + tumor.RR)]
inputDT[is.na(tumor.SR.ratio), tumor.SR.ratio := 0]

inputDT[, control.coverage := control.DR + control.DV + control.RR + control.RV]
inputDT[, tumor.coverage := tumor.DR + tumor.DV + tumor.RR + tumor.RV]

inputDT[, control.varFreq := (control.DV + control.RV) / control.coverage]
inputDT[is.na(control.varFreq), control.varFreq := 0]

inputDT[, tumor.varFreq := (tumor.DV + tumor.RV) / tumor.coverage]
inputDT[is.na(tumor.varFreq), tumor.varFreq := 0]

inputDT[donChr == accChr, distances := accPos - donPos]

# FILTERING
if(FilterInterChromosomalRearrangements){
  print(paste("Filtering only intrachromosomal rearrangements"))
  inputDT <- inputDT[donChr==accChr]
} 
filterDT <- inputDT[control.varFreq == 0 & tumor.varFreq >= TumorVarFreq & distances >= maxDist]
  
# stronger filtering for IMPRECISE variants
# Distance set because of huge FragmentSizes
# Control coverage set as filter for repeat regions , changed to 40 x
filterDT[IMPRECISE == T & distances > 1000000 & control.coverage <= CoverageFilter & MAPQ >= 30, remove := "yes"]
filterDT <- filterDT[is.na(remove)]

# REFORMAT
setorder(filterDT, donChr, donPos)
filterDT[PRECISE == T, Precision := "PRECISE"]
filterDT[IMPRECISE == T, Precision := "IMPRECISE"]
filterDT[, remove := NULL]
filterDT[, PRECISE := NULL]
filterDT[, IMPRECISE := NULL]
setcolorder(filterDT, c("donChr", "donPos", "accChr", "accPos", "PE", "SR", "Precision"))
  
print("Writing Output")
write.table(filterDT,paste(name,"/results/Delly/",name,".breakpoints.filtered.tab",sep=""),quote=F,row.names=F,col.names=T,sep='\t')


