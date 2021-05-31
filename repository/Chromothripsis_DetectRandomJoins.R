#!/usr/bin/Rscript

##########################################################################################
##
## Chromothripsis_DetectRandomJoins.R
##
## Check the random order and orientation of chromothriptic fragments.
##
##########################################################################################

message("\n###Random Fragment Order and Orientation###")
options(warn=-1)
suppressPackageStartupMessages(library(optparse))
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

# extract DELLY breakpoint list from rearrangement file
delly_breaks <- data.frame(input_list$donPos,input_list$accPos,as.character(input_list$Type),stringsAsFactors=F)
colnames(delly_breaks) <- c('Start','End','Type')
breakpoint_list <- sort(c(delly_breaks$Start,delly_breaks$End))

# Randomness of Fragment Joins
typeTab=c(sum(delly_breaks$Type=="3to5"),sum(delly_breaks$Type=="5to3"),sum(delly_breaks$Type=="3to3"),sum(delly_breaks$Type=="5to5"))
test <- chisq.test(typeTab,p=c(0.25,0.25,0.25,0.25))
#message(paste0("Randomness of Fragment Joins: p = ",signif(test$p.value,2)))

# Randomness of Fragment Order
dist <- c()
for(i in 1:nrow(delly_breaks)){
  b1 <- which(breakpoint_list == delly_breaks$Start[i])
  b2 <- which(breakpoint_list == delly_breaks$End[i])
  dist <- c(dist,abs(b2-b1))
}
sample_mean <- mean(dist)

means <- c()
for(i in 1:1000){
  random_breaks <- sample(seq(1,length(breakpoint_list)))
  rand.dist <- c()
  for(j in seq(2,length(breakpoint_list),by=2)){
    rand.dist <- c(rand.dist,abs(random_breaks[j]-random_breaks[j-1]))
  }
  means <- c(means,mean(rand.dist))
}

tmp = hist(means,freq=F,plot=F)
ymax = 1.1*max(tmp$density)
xmax = round(1.1*max(tmp$breaks))

# generate plots
if(opt$format=="tif") { tiff(paste0(opt$name,"/results/Chromothripsis/Chr",opt$chrom,"/",opt$name,".chr",opt$chr,".RandomFragmentOrder.tif"),1600,1600,res=200)
} else emf(paste0(opt$name,"/results/Chromothripsis/Chr",opt$chrom,"/",opt$name,".chr",opt$chr,".RandomFragmentOrder.emf"),bg="white",width=9,height=9,coordDPI = 200)
	par(oma=c(3,3,2,2),axes=F)

	plot(seq(0,xmax,by=0.5),dnorm(x=seq(0,xmax,by=0.5),mean=mean(means),sd=sd(means)),xlim=c(0,xmax),ylim=c(-0.1,ymax+0.1),
	     axes=F,type='l',las=1,frame.plot=F,xlab='',ylab='',col='red',lwd=2)
	hist(means,freq=F,add=T,lty=2)
	segments(sample_mean,0,sample_mean,ymax,col='blue',lwd=2)
	p_value <- 2*min(pnorm(sample_mean,mean=mean(means),sd=sd(means)),1-pnorm(sample_mean,mean=mean(means),sd=sd(means)))

	axis(side=2,las=1,at=seq(0,ymax+0.05,by=0.05),labels=seq(0,ymax+0.05,by=0.05),line=0,lwd=2,cex.axis=2)
	mtext(side=2,las=3,at=0.5*ymax,line=5,c('Density'),cex=2.2)

	text(x=sample_mean-5,y=0.9*ymax,labels=paste0('p = ',signif(p_value,2)),cex=2)

	segments(c(0,0,xmax),c(-0.01,-0.01,-0.01),c(xmax,0,xmax),c(-0.01,-0.02,-0.02),cex=2.2)
	arrows(x0=c(0.5*xmax,0.5*xmax),y0=c(-0.03,-0.03),x1=c(0,xmax),y1=c(-0.03,-0.03),length=c(0.1,0.1),lwd=2)
	text(x=c(0.1*xmax,0.9*xmax),y=c(-0.05,-0.05),labels=c('ordered','random'),cex=2.2)
garbage <- dev.off()


if(opt$format=="tif") { tiff(paste0(opt$name,"/results/Chromothripsis/Chr",opt$chrom,"/",opt$name,".chr",opt$chr,".RandomFragmentJoins.tif"),width=1600,height=700,res=100)
} else emf(paste0(opt$name,"/results/Chromothripsis/Chr",opt$chrom,"/",opt$name,".chr",opt$chr,".RandomFragmentJoins.emf"),bg="white", width=16,height=9,coordDPI = 200)
  pie(typeTab,labels="",col=c("#76323F","#57BC90","#D9B310","#438496"),border="white")
  legend(1,0.4,legend=c("3'-to-5' deletion-type","5'-to-3' duplication-type","3'-to-3' inversion-type",
                          "5'-to-5' inversion-type"),fill=c("#76323F","#57BC90","#D9B310","#438496"),bty="n",cex=2)
  text(-1.8,0,paste0("Goodness-of-fit test (n = ",nrow(delly_breaks),")\np = ",round(test$p.value,2)),cex=2)
garbage <- dev.off()