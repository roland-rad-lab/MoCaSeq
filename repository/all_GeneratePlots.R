#!/usr/bin/Rscript

##########################################################################################
##
## all_GeneratePlots.R
##
## Plots datapoints and segments for LOH and CNV
##
##########################################################################################

library(GenomicRanges)
library(devEMF)

DefineChromSizes = function(GenomeVersion)
{
  if(species=="Human")
  {
  chrom.sizes = c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,
                  138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,
                  83257441,80373285,58617616,64444167,46709983,50818468)
  names(chrom.sizes) = c(1:22)
  }

  if(species=="Mouse")
  {
  chrom.sizes = c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,
                  130694993,122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566)
  names(chrom.sizes) = c(1:19)
  }
  return(chrom.sizes)
}

FindCordsChromNames = function(ChromBorders)
{
ChromNamePos = c()
for ( i in 1:length(ChromBorders) )
{
  if(i==1)
  {
    ChromNamePos = ChromBorders[i]/2
  } else {
    ChromNamePos = c(ChromNamePos,(ChromBorders[i-1] + ChromBorders[i])/2)
  }
}
names(ChromNamePos) = names(ChromBorders)
return(ChromNamePos)
}


ConvertGenomicCords = function(dat,chrom.sizes,Start,CopyNumber,Chromosome)
{
  Outlist = list()
  cn <- data.frame()
  borders <- c()
  ChromBorders = c()
  FirstPosition = c()
  last = 0
  for(i in names(chrom.sizes))
  {
    cur <- dat[dat[,Chromosome]==i,c(Start,CopyNumber)]
    if (cur > 0) 
    {
       cn <- rbind(cn,data.frame(Chromosome = i,position=cur[,Start]+last,copy=cur[,CopyNumber]))
    } else {
       cn <- rbind(cn,data.frame(Chromosome = i,position=cur[,Start]+last,copy=cur[,CopyNumber]))
    }
    borders <- c(borders,last)
    last = last + chrom.sizes[i]
    ChromBorders = c(ChromBorders,last)
    FirstPosition = c(FirstPosition,min(dat[dat[,Chromosome]==i,Start]))
  }
  names(ChromBorders) = names(chrom.sizes)
  Outlist[["CN"]]=cn
  Outlist[["ChromosomeBorders"]] = ChromBorders
  Outlist[["FirstChromPosition"]] = FirstPosition
  Outlist[["NonProcessed"]] = dat
  return(Outlist)
}


SetVariableNames = function(method)
{
if(method=="Copywriter")
  {
    Start <<- "Start"
    CopyNumber <<- "log2Ratio"
    Chromosome <<- "Chrom"
  }
if(method=="HMMCopy")
  {
    Start <<- "Start"
    CopyNumber <<- "log2Ratio"
    Chromosome <<- "Chrom"
  }
if(method=="LOH")
  {
    Start <<- "Pos"
    CopyNumber <<- "Plot_Freq"
    Chromosome <<- "Chrom"
  }
if(method=="CNVKit")
  {
    Start <<- "start"
    CopyNumber <<- "log2"
    Chromosome <<- "chromosome"
  }
}


SetVariableNamesSegments = function(method)
{
if(method=="Copywriter")
  {
  SegmentStart <<- "Start"
  SegmentStop <<- "End"
  SegmentMean <<- "Mean"
  SegmentChromosome <<- "Chrom"
  }
if(method=="HMMCopy")
  {
  SegmentStart <<- "Start"
  SegmentStop <<- "End"
  SegmentMean <<- "Mean"
  SegmentChromosome <<- "Chrom"
if(method=="cnvkit")
  {
  SegmentStart <<- "start"
  SegmentStop <<- "end"
  SegmentMean <<- "log2"
  SegmentChromosome <<- "chromosome"
  }
}


ProcessCountData = function(countdata="",chrom.sizes=chrom.sizes,method="")
{
  Outlist = list()
  countdata = read.table(countdata,header=T,sep="\t")
  SetVariableNames(method)
  Outlist = ConvertGenomicCords(countdata,chrom.sizes,Start,CopyNumber,Chromosome)
  return(Outlist)
}


ProcessSegmentData = function(segmentdata="",chrom.sizes=chrom.sizes,method="")
{
  Outlist = list()
  FirstPosition = c()
  segmentdata = read.table(segmentdata,header=T,sep="\t")
  SetVariableNamesSegments(method)
  cnSeg <- data.frame()
  borders <- c()
  last = 0
  for(i in names(chrom.sizes))
  {
    cur <- segmentdata[segmentdata[,SegmentChromosome]==i,c(SegmentChromosome,SegmentStart,SegmentStop,SegmentMean)]
    cur[cur[,SegmentMean]<=(-5),SegmentMean]=-4.9
    cur[cur[,SegmentMean]>=(5),SegmentMean]=4.9
    cnSeg <- rbind(cnSeg,data.frame(Chromosome=cur[,SegmentChromosome],start=cur[,SegmentStart]+last,
                    stop=cur[,SegmentStop]+last,copy=cur[,SegmentMean]))
    borders <- c(borders,last)
    last = last + chrom.sizes[i]
  }
  Outlist[["CN"]] = cnSeg
  Outlist[["NonProcessed"]] = segmentdata
  return(Outlist)
}


SetYAxis = function(y_axis)
{
  if(y_axis=="CNV_5")
  {
    ylim <<- c(-5.5,5.5)
    ypos <<- seq(-5,5,1)
    ylabels <<- seq(-5,5,1)
    y_output <<- ".5."
    yborder <<- 5
  }
  if(y_axis=="CNV_2")
  {
    ylim <<- c(-2.5,2.5)
    ypos <<- seq(-2,2,1)
    ylabels <<- seq(-2,2,1)
    y_output <<- ".2."
    yborder <<- 2
  }
  if(y_axis=="LOH")
   {
     ylim <<- c(-0.25,1.25)
     ypos <<- c(0,0.5,1)
     ylabels <<- c(0,0.5,1)
     y_output <<- "."
   }
}


SetXAxis = function(chrom.sizes,chromosome,SliceStart,SliceStop,Organism)
{
ValueRange = seq(0,max(chrom.sizes),20000000)
ValueRange = c(ValueRange,max(ValueRange)+20000000)
if((SliceStart=="") & (SliceStop==""))
  {
    Xmax = as.numeric(chrom.sizes[chromosome])
    Xmax <<- min(ValueRange[ValueRange>=Xmax])
    Xmin <<- 0
    StepSize <<- 20000000
  } else {
    Xmax <<- SliceStop
    Xmin <<- SliceStart
    StepSize <<- 0.1* (Xmax-Xmin)
  }
}


plotGlobalRatioProfile = function(cn=cn,ChromBorders=ChromBorders,cnSeg="",samplename="",method="",toolname="",normalization="",y_axis="",Transparency=20, Cex=0.2,outformat="")
{
  SetYAxis(y_axis)
  SetVariableNamesSegments(toolname)
  SetVariableNames(toolname)
  ChromNamePos = FindCordsChromNames(ChromBorders)
  ChromBorders=c(0,ChromBorders)
  YaxisPosition = min(ChromBorders)-130000000
  LastEntry = length(ChromBorders)
  if(normalization != "")
  {
    normalization=paste(".",normalization,sep="")
  }
  if(outformat=="emf")
  {
    emf(paste(samplename,".",method,".",toolname,normalization,y_output,"emf",sep=""), bg="white", width=25,height=10)
  }
  if(outformat=="pdf")
  {
    pdf(paste(samplename,".",method,".",toolname,normalization,y_output,"pdf",sep=""),width=25,height=10)
  }
  par(mar=c(4, 4, 0, 0))
  plot(cn$position,cn$copy,pch=20,cex=Cex,
      ylim=ylim,xaxt="n",yaxt="n",bty="n",col=paste("#000000",Transparency,sep=""),yaxs="i",ylab="",xlab="",yaxt="n")
  segments(min(ChromBorders),0,ChromBorders[LastEntry],col="#A9A9A9",lwd=2)
  axis(2,las=1,pos=YaxisPosition, outer=T, at=ypos,labels=ylabels)
  axis(1,las=1, labels=rep("",length(ChromBorders)),at=ChromBorders)
  mtext(side=1,line=1,at=ChromNamePos,names(ChromNamePos))
  if(method=="CNV")
  {
    segments(cnSeg$start,cnSeg$copy,cnSeg$stop,cnSeg$copy,col="#FF4D4D",lwd=4)
  }
  ChromBordersReduced = ChromBorders[-c(1,LastEntry)]
  if(method == "CNV")
  {
    segments(ChromBordersReduced,-yborder,ChromBordersReduced,yborder,lty=3,col="grey40",lwd=0.9)
  }
  if(method =="LOH")
  {
    segments(ChromBordersReduced,0,ChromBordersReduced,1,lty=3,col="grey40",lwd=0.9)
  }
  dev.off()
}


plotChromosomalRatioProfile = function(cn=cn,chrom.sizes,cnSeg="",samplename="",chromosome="",method="",toolname="",normalization="",y_axis="",SliceStart="",SliceStop="",Transparency=70, Cex=0.3,outformat="")
{
  SetYAxis(y_axis)
  SetVariableNamesSegments(toolname)
  SetVariableNames(toolname)
  cn = cn[cn[,Chromosome]==chromosome,]
  SetXAxis(chrom.sizes,chromosome,SliceStart,SliceStop)
  if((SliceStart!="") & (SliceStop!=""))
  {
    cn = cn[cn[,Start]>=SliceStart & cn[,Start]<=SliceStop,]
  }
  if((SliceStart=="") & (SliceStop==""))
  {
    SliceStart=0
    SliceStop=chrom.sizes[chromosome]
  }
  if(normalization != "")
  {
    normalization=paste(".",normalization,sep="")
  }
  if(outformat == "emf")
  {
    emf(paste(samplename,"_Chromosomes/",samplename,".Chr",chromosome,".",method,".",toolname,normalization,y_output,"emf",sep=""), bg="white", width=25,height=10)
  }
  if(outformat == "pdf")
  {
    pdf(paste(samplename,"_Chromosomes/",samplename,".Chr",chromosome,".",method,".",toolname,normalization,y_output,"pdf",sep=""),width=25,height=10)
  }

  par(mar=c(4, 4, 0, 0))
  plot(cn[,Start],cn[,CopyNumber],pch=20,cex=Cex,xlim=c(Xmin,Xmax),ylim=ylim,xaxt="n",yaxt="n",bty="n",col=paste("#000000",Transparency,sep=""),ylab="",xlab="",yaxt="n")
  if(method=="CNV")
  {
    cnSeg = cnSeg[cnSeg[,SegmentChromosome]==chromosome,]
    cnSeg[cnSeg[,SegmentMean]<=(-5),SegmentMean]=-4.9
    cnSeg[cnSeg[,SegmentMean]>=(5),SegmentMean]=4.9
    segments(cnSeg[1,SegmentStart],0,cnSeg[nrow(cnSeg),SegmentStop],col="#5C5C5C",lwd=4)
    segments(cnSeg[,SegmentStart],cnSeg[,SegmentMean],
             cnSeg[,SegmentStop],cnSeg[,SegmentMean],col="#FFFFFF",lwd=6)
    usr = par("usr")
    clip(SliceStart,SliceStop, usr[3], usr[4])
    segments(cnSeg[,SegmentStart],cnSeg[,SegmentMean],
             cnSeg[,SegmentStop],cnSeg[,SegmentMean],col="#FF4D4D",lwd=6)
  }
  YAxisScaling = (Xmax-Xmin)*0.015
  axis(1,at=seq(Xmin,Xmax,StepSize),labels=(seq(Xmin,Xmax,StepSize) / 1000000))
  axis(2,las=1,pos=Xmin-YAxisScaling, outer=T, at=ypos,labels=ylabels)
  title(main = paste(chromosome,sep=""), line = -3, font.main = 1)
  dev.off()
}