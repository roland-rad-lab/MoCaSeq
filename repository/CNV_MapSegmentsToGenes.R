#!/usr/bin/Rscript

##########################################################################################
##
## CNV_MapSegmentsToGenes.R
##
## Takes the segment file from either HMMCopy or Copywriter and maps them to genes.
##
##########################################################################################

args <- commandArgs(TRUE)

name = args[1]
species = args[2]
genecode_file = args[3]
method = args[4]
resolution = args[5]
CGC = args[6]
TruSight = args[7]

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))

genesDT = readRDS(genecode_file)

genesGR <- makeGRangesFromDataFrame(genesDT, keep.extra.columns = T)

AnnotateSegment <- function(segDF)
{
  segDF[,"chr"]=paste0("chr",segDF[,"chr"])
  segGR <- makeGRangesFromDataFrame(segDF, keep.extra.columns = T)
  hits <- findOverlaps(segGR, genesGR)
  
  #returnDat <- genesDT[subjectHits(hits), .(chr, start, end, geneID)] # only geneID
  returnDat <- genesDT[subjectHits(hits), .(chr, start, end, geneName)] # more stuff
  return(data.frame(returnDat))
}

if (method=="Copywriter")
{
	segment=paste(name,"/results/",method,"/",name,".",method,".segments.Mode.txt",sep="")
} else if (method=="HMMCopy") {
	segment=paste(name,"/results/",method,"/",name,".",method,".",resolution,".segments.txt",sep="")
}

cnv=data.frame(Name=NULL,Chrom=NULL, Start=NULL, End=NULL, Mean=NULL,Gene=NULL)
segment=read.delim(segment)

for (i in 1:nrow(segment))
{
	temp=NULL
	temp=as.data.frame(segment[i,c("Chrom","Start","End")])
	colnames(temp)=c("chr","start","end")
	results=AnnotateSegment(temp)
	if (nrow(results) > 0) 
		{
		colnames(results)=c("Chrom", "Start", "End", "Gene")
		results[,"Chrom"]=gsub("chr","",results[,"Chrom"])
		results$Mean=segment[i,"Mean"]
		results$Name=name
		results=results[,c("Name", "Chrom", "Start", "End", "Mean","Gene")]
		cnv=rbind(cnv,results)
		}
}

if (species=="Mouse")
{
	segment=makeGRangesFromDataFrame(segment,keep.extra.columns=T)
	ncruc.gr=GRanges(4, IRanges(89311040, 89511040))
	olaps=findOverlaps(segment,ncruc.gr)
	ncruc=data.frame(pintersect(segment[queryHits(olaps)], ncruc.gr[subjectHits(olaps)]))[,c("seqnames","start","end","Mean")]
	ncruc$Gene="Cdkn2_ncruc"
	ncruc$Name=name
	colnames(ncruc)=c("Chrom", "Start", "End", "Mean", "Gene", "Name")
	ncruc=ncruc[,c("Name","Chrom", "Start", "End", "Mean", "Gene")]
	cnv=rbind(cnv,ncruc)
}

cnv=tbl_df(cnv) %>% mutate_each(as.character)

#print("For all genes on segment border, pick the absolute largest mean.")

cnv = cnv %>% 
group_by(Gene) %>% 
arrange(desc(abs(as.numeric(Mean))),.by_group=T) %>% 
filter(row_number()==1) %>% 
arrange(as.numeric(Chrom),as.numeric(Start),as.numeric(End))

if (method=="Copywriter")
{
	write.table(cnv,paste(name,"/results/",method,"/",name,".",method,".genes.Mode.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
} else if (method=="HMMCopy") {
	write.table(cnv,paste(name,"/results/",method,"/",name,".",method,".",resolution,".genes.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
}

CGC=read.delim(CGC,header=T,sep="\t")

cnv_cgc = cnv %>%
filter(Gene %in% CGC[,1])

if (method=="Copywriter")
{
	write.table(cnv_cgc,paste(name,"/results/",method,"/",name,".",method,".genes.Mode.CGC.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
} else if (method=="HMMCopy") {
	write.table(cnv_cgc,paste(name,"/results/",method,"/",name,".",method,".",resolution,".genes.CGC.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
}

cnv_cgc = cnv_cgc %>%
filter(abs(as.numeric(Mean)) > 0.75)

if (method=="Copywriter")
{
	write.table(cnv_cgc,paste(name,"/results/",method,"/",name,".",method,".genes.Mode.OnlyImpact.CGC.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
} else if (method=="HMMCopy") {
	write.table(cnv_cgc,paste(name,"/results/",method,"/",name,".",method,".",resolution,".genes.OnlyImpact.CGC.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
}

if (species=="Human")
{
	TruSight=read.delim(TruSight,header=T,sep="\t")

	cnv_ts = cnv %>%
	filter(Gene %in% TruSight[,1])

	if (method=="Copywriter")
	{
		write.table(cnv_ts,paste(name,"/results/",method,"/",name,".",method,".genes.Mode.TruSight.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
	} else if (method=="HMMCopy") {
		write.table(cnv_cgc,paste(name,"/results/",method,"/",name,".",method,".",resolution,".genes.TruSight.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
	}

	cnv_ts = cnv_ts %>%
	filter(abs(as.numeric(Mean)) > 0.75)

	if (method=="Copywriter")
	{
		write.table(cnv_ts,paste(name,"/results/",method,"/",name,".",method,".genes.Mode.OnlyImpact.TruSight.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
	} else if (method=="HMMCopy") {
		write.table(cnv_ts,paste(name,"/results/",method,"/",name,".",method,".",resolution,".genes.OnlyImpact.TruSight.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
	}
}