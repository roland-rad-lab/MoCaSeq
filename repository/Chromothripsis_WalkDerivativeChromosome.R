#!/usr/bin/Rscript

##########################################################################################
##
## Chromothripsis_WalkDerivateChromosome.R
##
## Analyse the alternating behavior of 5' and 3' joins.
##
##########################################################################################

message("\n###HT Sequencer###")
options(warn=-1)
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(randtests))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(devEMF))

option_list = list(
  make_option(c("-i", "--input"),type="character",default=NULL,help="rearrangement list"),
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
if(opt$chrom %in% input_list$donChr){
  input_list = input_list[input_list$donChr==opt$chrom & input_list$accChr==opt$chrom,]
} else{
  stop("The specified chromosome couldn't be found in the rearrangement list.",call.=FALSE)
}
if(nrow(input_list)<4) stop("There are too few (< 4) rearrangements on the specified chromosome. No analysis possible.",call.=FALSE)
if(!opt$format %in% c("tif","emf")) stop("The output format has to be either tif or emf.",call.=FALSE)


# main program

# DELLY Breakpoint List
delly_breaks <- data.frame(input_list$donPos,input_list$accPos,as.character(input_list$Type),stringsAsFactors=F)
colnames(delly_breaks) <- c('Start','End','Type')

# Walk the Derivative Chromosome
starts <- data.frame(delly_breaks$Start,substr(delly_breaks$Type,1,1),stringsAsFactors=F)
colnames(starts) <- c('bp','Type')
ends <- data.frame(delly_breaks$End,substr(delly_breaks$Type,4,4),stringsAsFactors=F)
colnames(ends) <- c('bp','Type')

seq53 <- rbind(starts,ends)
seq53 <- seq53[order(seq53$bp),]
seq53$Type = as.numeric(seq53$Type)

test <- runs.test(seq53$Type,threshold=4,alternative='r')
#print(test)


# generate plot
boxes <- data.frame(x=c(1:nrow(seq53)),color=seq53$Type)

del <- data.frame()
for(i in 1:nrow(input_list[input_list$Type == '3to5',])){
  del <- rbind(del,data.frame(donPos = which(seq53[,1] == input_list[input_list$Type == '3to5',2][i]), accPos=which(seq53[,1] == input_list[input_list$Type == '3to5',4][i])))
}
dup <- data.frame()
for(i in 1:nrow(input_list[input_list$Type == '5to3',])){
  dup <- rbind(dup,data.frame(donPos = which(seq53[,1] == input_list[input_list$Type == '5to3',2][i]), accPos=which(seq53[,1] == input_list[input_list$Type == '5to3',4][i])))
}
inv1 <- data.frame()
for(i in 1:nrow(input_list[input_list$Type == '3to3',])){
  inv1 <- rbind(inv1,data.frame(donPos = which(seq53[,1] == input_list[input_list$Type == '3to3',2][i]), accPos=which(seq53[,1] == input_list[input_list$Type == '3to3',4][i])))
}
inv2 <- data.frame()
for(i in 1:nrow(input_list[input_list$Type == '5to5',])){
  inv2 <- rbind(inv2,data.frame(donPos = which(seq53[,1] == input_list[input_list$Type == '5to5',2][i]), accPos=which(seq53[,1] == input_list[input_list$Type == '5to5',4][i])))
}

if(test$p.value < 0.001){
  pText = paste0("p < 10E",ceiling(log10(test$p.value)))
} else{
  pText = paste0("p = ",signif(test$p.value,2))
}

p <- ggplot(boxes,aes()) + ylim(-2,2.5) + geom_rect(aes(xmin = x-1,xmax = x,ymin = 0,ymax = 0.5, fill = as.factor(color)),color='grey') + scale_fill_manual(values=c('#DE4E4E','#D9D9D9'),labels = c("3' End","5' End"),guide = guide_legend(title = "Segment End")) #Define box colors at scale_fill_manual
  if(nrow(del)!=0) p <- p + geom_segment(aes(x = donPos-0.5, y = 1.5, xend = donPos-0.5, yend = 0.5, colour = "Deletion-type"), data = del,size=0.6) +
                            geom_segment(aes(x = accPos-0.5, y = 1.5, xend = accPos-0.5, yend = 0.5, colour = "Deletion-type"), data = del,size=0.6) +
                            geom_curve(aes(x = donPos-0.5, y = 1.5, xend = accPos-0.5, yend = 1.5, colour = "Deletion-type"), data = del,curvature=-0.1,size=0.6)

  if(nrow(dup)!=0) p <- p + geom_segment(aes(x = donPos-0.5, y = 1.5, xend = donPos-0.5, yend = 0.5, colour = "Duplication-type"), data = dup,size=0.6) +
                            geom_segment(aes(x = accPos-0.5, y = 1.5, xend = accPos-0.5, yend = 0.5, colour = "Duplication-type"), data = dup,size=0.6) +
                            geom_curve(aes(x = donPos-0.5, y = 1.5, xend = accPos-0.5, yend = 1.5, colour = "Duplication-type"), data = dup,curvature=0.1,size=0.6)

  if(nrow(inv1)!=0) p <- p + geom_segment(aes(x = donPos-0.5, y = -1, xend = donPos-0.5, yend = 0, colour = "Inversion-type 1"), data = inv1,size=0.6) +
                             geom_segment(aes(x = accPos-0.5, y = -1, xend = accPos-0.5, yend = 0, colour = "Inversion-type 1"), data = inv1,size=0.6) +
                             geom_curve(aes(x = donPos-0.5, y = -1, xend = accPos-0.5, yend = -1, colour = "Inversion-type 1"), data = inv1,curvature=-0.1,size=0.6)

  if(nrow(inv2)!=0) p <- p + geom_segment(aes(x = donPos-0.5, y = -1, xend = donPos-0.5, yend = 0, colour = "Inversion-type 2"), data = inv2,size=0.6) +
                             geom_segment(aes(x = accPos-0.5, y = -1, xend = accPos-0.5, yend = 0, colour = "Inversion-type 2"), data = inv2,size=0.6) +
                             geom_curve(aes(x = donPos-0.5, y = -1, xend = accPos-0.5, yend = -1, colour = "Inversion-type 2"), data = inv2,curvature=0.1,size=0.6)

  p <- p + scale_color_manual(values=c("#76323F","#57BC90","#D9B310","#438496"),guide = guide_legend(title = "Rearrangement Type"))

  p <- p + geom_segment(aes(x = 0, y = 1.5, xend = nrow(seq53), yend = 1.5), color = 'black',size=1.2) +
           geom_segment(aes(x = 0, y = -1, xend = nrow(seq53), yend = -1), color = 'black',size=1.2) +
           annotate("text", x = c(-5,-5,-5), y = c(1.5,-1,0.25), label = c("Non-inverted\norientation","Inverted\norientation",pText),size=4)

  p <- p + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),
                 axis.title.y=element_blank(),panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),plot.background=element_blank())

if(opt$format=="tif"){ tiff(paste0(opt$name,"/results/Chromothripsis/Chr",opt$chrom,"/",opt$name,".chr",opt$chrom,".53Seq.tif"),width=2600,height=900,res=200)
} else emf(paste0(opt$name,"/results/Chromothripsis/Chr",opt$chrom,"/",opt$name,".chr",opt$chrom,".53Seq.emf"),bg="white", width=12,height=6)
  plot(p)
garbage <- dev.off()