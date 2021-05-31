#!/usr/bin/Rscript

##########################################################################################
##
## Chromothripsis_PlotRearrangementGraph.R
##
## Produces a rearrangement graph.
##
##########################################################################################

message("\n###Rearrangement Graph###")
options(warn=-1)
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-i", "--input"),type="character",default=NULL,help="rearrangement list"),
  make_option(c("-d", "--copynumber"),type="character",default=NULL,help="file containing copy number datapoints"),
  make_option(c("-c", "--chrom"),type="character",default=NULL,help="chromosome to analyze"),
  make_option(c("-n", "--name"),type="character",default="Sample1",help="sample name [default = %default]"),
  make_option(c("-f", "--format"),type="character",default="tif",help="output format (tif|emf) [default = %default]")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$input)){
  print_help(opt_parser)
  stop("You have to specify a rearrangement list.",call.=FALSE)
} else{
  tryCatch({
    input_list <- read.csv(opt$input,header=T,sep="\t",stringsAsFactors=F)
  },
  error = function(e){
    stop("Error reading rearrangement list.")
  })
}
if(is.null(opt$copynumber)){
  print_help(opt_parser)
  stop("You have to specify a copy number data file.", call.=FALSE)
} else{
  tryCatch({
    CNdatapoints <- read.csv(opt$copynumber,header=T,sep="\t",stringsAsFactors=F)
  },
  error = function(e){
    stop("Error reading copy number data file.")
  })
}
if(opt$chrom %in% input_list$donChr){
  input_list = input_list[input_list$donChr==opt$chrom & input_list$accChr==opt$chrom,]
} else{
  stop("The specified chromosome couldn't be found in the rearrangement list.",call.=FALSE)
}
if(nrow(input_list)<4) stop("There are too few (< 4) rearrangements on the specified chromosome. No plotting.",call.=FALSE)
if(!opt$format %in% c("tif","emf")) stop("The output format has to be either tif or emf.",call.=FALSE)

# main program
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("devEMF")) install.packages("devEMF")
suppressMessages(library(ggplot2))
suppressMessages(library(devEMF))

table_draw <- CNdatapoints[CNdatapoints$Chrom==opt$chrom & !is.na(CNdatapoints$log2Ratio),]
Max = max(table_draw$End)

# DELLY Breakpoint List
delly_breaks <- data.frame(input_list$donPos,input_list$accPos,as.character(input_list$Type),stringsAsFactors=F)
colnames(delly_breaks) <- c('Start','End','Type')

delly_breaks$x = 0
delly_breaks$y = 0
delly_breaks$xend = 0
delly_breaks$yend = 0
delly_breaks$Curvature=0
delly_breaks$Rearrangement="None"

Positions = c(delly_breaks[1:nrow(delly_breaks),c("Start")],delly_breaks[1:nrow(delly_breaks),"End"])
SVType = c(delly_breaks[1:nrow(delly_breaks),"Type"],delly_breaks[1:nrow(delly_breaks),"Type"])
delly_formated = data.frame(Position=Positions,SVType=SVType)

delly_formated$Rearrangement = "None"

delly_formated[delly_formated$SVType=="3to5","Rearrangement"]="Deletion-type"
delly_formated[delly_formated$SVType=="5to3","Rearrangement"]="Duplication-type"
delly_formated[delly_formated$SVType=="3to3","Rearrangement"]="Inversion-type 1"
delly_formated[delly_formated$SVType=="5to5","Rearrangement"]="Inversion-type 2"

delly_formated$y=0
delly_formated$yend=0
delly_formated[delly_formated$Rearrangement %in% c("Deletion-type","Duplication-type"),"yend"] = 2
delly_formated[delly_formated$Rearrangement %in% c("Inversion-type 1","Inversion-type 2"),"yend"] = -2

delly_breaks[delly_breaks$Type=="3to5","Rearrangement"]="Deletion-type"
delly_breaks[delly_breaks$Type=="5to3","Rearrangement"]="Duplication-type"
delly_breaks[delly_breaks$Type=="3to3","Rearrangement"]="Inversion-type 1"
delly_breaks[delly_breaks$Type=="5to5","Rearrangement"]="Inversion-type 2"

delly_breaks[delly_breaks$Type %in% c("3to5","5to3"),"y"] = 2
delly_breaks[delly_breaks$Type %in% c("3to3","5to5"),"y"] = -2
delly_breaks[delly_breaks$Type %in% c("3to5","5to3"),"Curvature"] = -0.1
delly_breaks[delly_breaks$Type %in% c("3to3","5to5"),"Curvature"] = 0.1

delly_breaks$Rearrangement = as.factor(delly_breaks$Rearrangement)
delly_formated$Rearrangement = as.factor(delly_formated$Rearrangement)

set1 = delly_breaks[delly_breaks$y==2,]
set2 = delly_breaks[delly_breaks$y==-2,]

set1$Type = as.factor(set1$Type)
set2$Type = as.factor(set2$Type)

# generate plot
var = ggplot(table_draw, aes(x=Start, y=log2Ratio)) + geom_point(size=0.5)
var = var + geom_segment(aes(x = Position , y = y, xend = Position , yend = yend ,colour=Rearrangement), data = delly_formated) +
  geom_curve(aes(x = Start, y = y, xend = End, yend = y,colour=Rearrangement), curvature = -0.1, data=set1) +
  geom_curve(aes(x = Start, y = y, xend = End, yend = y,colour=Rearrangement), curvature = 0.1, data=set2) +
  scale_color_manual(values=c("#76323F","#57BC90","#D9B310","#438496")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
  xlab("Chromosomal Position") + xlim(0,Max) +
  ylab("Log2 Ratio")+ylim(c(-4,4))

if(opt$format=="tif"){ tiff(paste0(opt$name,"/results/Chromothripsis/Chr",opt$chrom,"/",opt$name,".chr",opt$chrom,".RearrangementGraph.tif"),1000,500)
} else emf(paste0(opt$name,"/results/Chromothripsis/Chr",opt$chrom,"/",opt$name,".chr",opt$chrom,".RearrangementGraph.emf"),bg="white", width=25,height=25)
  plot(var)
garbage <- dev.off()