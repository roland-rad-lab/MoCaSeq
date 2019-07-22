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
method = args[3]
resolution = args[4]

suppressMessages(library(biomaRt))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(GenomicRanges))

if (method=="Copywriter")
{
	segment=paste(name,"/results/",method,"/",name,".",method,".segments.Mode.txt",sep="")
} else if (method=="HMMCopy") {
	segment=paste(name,"/results/",method,"/",name,".",method,".",resolution,".segments.txt",sep="")
}

if (species=="Mouse")
{
	dataset="mmusculus_gene_ensembl"
} else if (species=="Human") {
	dataset="hsapiens_gene_ensembl"
}

mart=useMart(biomart = 'ensembl', dataset = dataset)

cnv=data.frame(Name=NULL,Chrom=NULL, Start=NULL, End=NULL, Mean=NULL,Gene=NULL)
segment=read.delim(segment)

for (i in 1:nrow(segment))
{
	print(paste("Getting information from segment ",i,"/",nrow(segment),".",sep=""))
	results=getBM(attributes=c("chromosome_name","start_position","end_position","external_gene_name"), 
		filters=c("chromosomal_region"),
		values=paste(segment[i,"Chrom"],":",segment[i,"Start"],":",segment[i,"End"],sep=""),
		mart=mart)
	if (nrow(results) > 0) 
		{
		colnames(results)=c("Chrom", "Start", "End", "Gene")
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