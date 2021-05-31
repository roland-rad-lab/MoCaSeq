#!/usr/bin/Rscript

##########################################################################################
##
## Chromothripsis_DetectBreakpointClustering.R
##
## Compare the distribution of observed inter-breakpoint-distances to a theoretial exponential distribution.
##
##########################################################################################

message("\n##Breakpoint Cluster Analysis###")
options(warn=-1)
suppressPackageStartupMessages(library(optparse))

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
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(devEMF))

# extract DELLY breakpoint list from rearrangement file
delly_breaks <- data.frame(input_list$donPos,input_list$accPos,as.character(input_list$Type),stringsAsFactors=F)
colnames(delly_breaks) <- c('Start','End','Type')
breakpoint_list <- sort(c(delly_breaks$Start,delly_breaks$End))

# Breakpoint Cluster Analysis
dist <- c()
for(i in 2:length(breakpoint_list)){
  dist <- c(dist,breakpoint_list[i]-breakpoint_list[i-1])
}
lamda = 1/mean(dist) # ML estimator
dist_hist <- hist(dist,breaks=30,plot=F) # observed distribution

# generate theoretical distribution
breaks_cdf <- pexp(dist_hist$breaks,rate=lamda)
null.probs <- rollapply(breaks_cdf, 2, function(x) x[2]-x[1])
test <- chisq.test(dist_hist$counts, p=null.probs, rescale.p=TRUE, simulate.p.value=TRUE)

if(is.na(test$p.value)){
  test$p.value <- 1
}

# generate plot
if(test$p.value < 0.001){
  pText = paste0("p < 10E",ceiling(log10(test$p.value)))
} else{
  pText = paste0("p = ",signif(test$p.value,2))
}
  
if(opt$format=="tif") { tiff(paste0(opt$name,"/results/Chromothripsis/Chr",opt$chrom,"/",opt$name,".chr",opt$chrom,".BreakpointCluster.tif"),1600,1600,res=200)
} else emf(paste0(opt$name,"/results/Chromothripsis/Chr",opt$chrom,"/",opt$name,".chr",opt$chrom,".breakpoint_cluster.emf"),bg="white", width=9,height=9,coordDPI = 200)
  par(oma=c(4,6,0,0))
  #ylabel = expression(paste(10^0,10^1,10^2,10^3,10^4,10^5,10^6,10^7,10^8,sep=" "))
  #ylabel=bquote(.(ylabel))
  ylabel = c('1','10','100','1,000','10,000','100,000','1,000,000','10,000,000','100,000,000')
  stripchart(list(Observed=log10(dist+1),Expected=log10(rexp(length(dist),rate=lamda)+1)),pch=20,vertical=T,method='jitter',las=1,frame.plot=F,xlab='',ylab='',xaxt='n',yaxt='n',ylim=c(-2,9),at=c(1,1.3))
  axis(side=2,las=1,at=log10(c(1,10,100,1000,10000,100000,1000000,10000000,100000000)),labels=ylabel,lwd=2.2)
  mtext(side=2,las=3,at=4,line=6,c('Distance between \n adjacent breakpoints (bp)'),cex=2.2)
  axis(side=1,at=c(1,1.3),labels=F,line=-3,lwd=2)
  text(x=c(1,1.3),y=-1.9,labels=c('Observed','Expected'),srt=45,adj=1,xpd=TRUE,cex=2.2)
  segments(1,8,1,8.3)
  segments(1.3,8,1.3,8.3)
  segments(1,8.3,1.3,8.3)
  text(x=1.15,y=8.7,labels=paste(pText,", (n=",length(dist),")",sep=""))
garbage <- dev.off()