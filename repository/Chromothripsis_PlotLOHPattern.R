#!/usr/bin/Rscript

##########################################################################################
##
## Chromothripsis_PlotLOHPattern.R
##
## Visualisation of the overlap between heterozygously deleted copy number segments and regions of LOH.
##
##########################################################################################

message("\n###LOH vs CN Comparison###")
options(warn=-1)
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-s", "--segments"),type="character",default=NULL,help="file containing copy number segments"),
  make_option(c("-d", "--copynumber"),type="character",default=NULL,help="file containing copy number datapoints"),
  make_option(c("-v", "--variants"),type="character",default=NULL,help="file containing variant positions"),
  make_option(c("-o", "--organism"),type="character",default=NULL,help="human|mouse"),
  make_option(c("-c", "--chrom"),type="character",default=NULL,help="chromosome to visualize"),
  make_option(c("-n", "--name"),type="character",default="Sample1",help="sample name [default = %default]"),
  make_option(c("-f", "--format"),type="character",default="tif",help="output format (tif|emf) [default = %default]")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

human_chroms <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,
                  135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,
                  50818468,156040895,57227415)
names(human_chroms)=c(1:22,"X","Y")
mouse_chroms <- c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,
                  122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566,171031299,91744698)
names(mouse_chroms)=c(1:19,"X","Y")


if(is.null(opt$segments)){
  print_help(opt_parser)
  stop("You have to specify a copy number segment file.", call.=FALSE)
} else{
  tryCatch({
    cnv_segments <- read.csv(opt$segments,header=T,sep="\t",stringsAsFactors=F)
  },
  error = function(e){
    stop("Error reading copy number segment file.")
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
if(is.null(opt$variants)){
  print_help(opt_parser)
  stop("You have to specify a variant file.", call.=FALSE)
} else{
  tryCatch({
    varFile <- read.csv(opt$variants,header=T,sep="\t",stringsAsFactors=F)
  },
  error = function(e){
    stop("Error reading variant file.")
  })
}
if(is.null(opt$organism)){
  print_help(opt_parser)
  stop("You have to specify an organism (either human or mouse).",call.=FALSE)
} else{
  if(is.null(opt$chrom)){
    print_help(opt_parser)
    stop("You have to specify a chromosome to visualize.",call.=FALSE)
  } else {
    if(opt$organism == "human" & opt$chrom %in% names(human_chroms)){
      chrom.sizes = human_chroms
      chr = opt$chrom
    }
    else if(opt$organism == "mouse" & opt$chrom %in% names(mouse_chroms)){
      chrom.sizes = mouse_chroms
      chr = opt$chrom
    }
    else{
      print_help(opt_parser)
      stop("The given organism and the specified chromsome do not match.", call.=FALSE)
    }
  }
}


# main program
if(!require("devEMF")) install.packages("devEMF")
suppressMessages(library(devEMF))

#print(varFile[1:4,])
table_draw <- CNdatapoints[CNdatapoints$Chrom==chr & !is.na(CNdatapoints$log2Ratio),]
break_draw <- cnv_segments[cnv_segments$Chrom==chr & abs(cnv_segments$Mean < -0.2),]
varFile = varFile[varFile$Chrom==chr,]
size = chrom.sizes[chr]

# assuming diploid copy number
ylim=c(-5,5)
#if(max(abs(range(table_draw$Mean))) <= 2) ylim=c(-2,2) else ylim=c(-5,5)

# sample datapoints for nice graphical representation
if(nrow(varFile)>nrow(table_draw)){ Ind = sample(1:nrow(varFile),nrow(table_draw))
} else Ind = 1:nrow(varFile)

if(opt$format=="tif"){ tiff(paste0(opt$name,"/results/Chromothripsis/Chr",opt$chrom,"/",opt$name,".chr",opt$chrom,".LOHvsCN.tif"),1000,500)
} else emf(paste0(opt$name,"/results/Chromothripsis/Chr",opt$chrom,"/",opt$name,".chr",opt$chrom,".LOHvsCN.emf"),bg="white",width=20,height=10)
	par(mfrow=c(2,1),mar=c(2,4,2,2))
    plot(table_draw$Start,table_draw$log2Ratio,pch=20,cex=0.4,col="#00000040",bty="n",xaxt="n",
         xlab="",ylab="Estimated CopyNumber",xlim=c(1,size),ylim=ylim,las=1)
    segments(cnv_segments[cnv_segments$Chrom==chr,'Start'],
         cnv_segments[cnv_segments$Chrom==chr,'Mean'],
         cnv_segments[cnv_segments$Chrom==chr,'End'],
         cnv_segments[cnv_segments$Chrom==chr,'Mean'],col='red')
    axis(side=1,c(1,seq(10e6,size,10e6)/1000000),at=c(1,seq(10e6,size,10e6)))

  par(mar=c(4,4,2,2))
	  plot(varFile$Pos[Ind],varFile$Plot_Freq[Ind],pch=20,cex=0.4,col="#00000040",bty="n",yaxt="n",xaxt="n",
	       ylab="B-allele frequency",xlim=c(1,size),xlab="Chromosome position (bp)")
    rect(break_draw$Start,min(varFile$Plot_Freq),break_draw$End,max(varFile$Plot_Freq),col='#66666610', border="white")
    axis(side=2,at=c(0,0.5,1),c(0,0.5,1),las=1)
    axis(side=1,c(1,seq(10e6,size,10e6)/1000000),at=c(1,seq(10e6,size,10e6)))

garbage <- dev.off()