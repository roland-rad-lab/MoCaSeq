#!/usr/bin/Rscript

##########################################################################################
##
## LOH_MapSegmentsToGenes.R
##
## Takes the segment file from Titan and maps it to genes.
##
##########################################################################################

args <- commandArgs(TRUE)

name = args[1]
species = args[2]
genecode_file = args[3]
CGC = args[4]
TruSight = args[5]

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))


genesDT = readRDS(genecode_file)

genesGR <- makeGRangesFromDataFrame(genesDT, keep.extra.columns = T)

AnnotateSegment <- function(segDF)
{
  segDF[,"chr"]=paste0("chr",segDF[,"chr"])
  segGR <- makeGRangesFromDataFrame(segDF, keep.extra.columns = T)
  hits <- findOverlaps(segGR, genesGR)
  
  #returnDat <- genesDT[subjectHits(hits), .(chr, start, end, geneID)] # only geneID
  returnDat <- genesDT[subjectHits(hits), .(chr, start, end, geneName, geneID)] # more stuff
  return(data.frame(returnDat))
}

segments=paste(name,"/results/Titan/run_ploidy2/",name,"_cluster01.segs.txt",sep="")
loh=data.frame(Name=NULL,Chrom=NULL, Start=NULL, End=NULL,TITAN=NULL,Gene=NULL)
segments=tbl_df(read.delim(segments))

segments = segments %>% filter(TITAN_call %in% c("ALOH","NLOH","DLOH","HOMD")) %>% filter(Length.snp. >= 10) %>% filter((End_Position.bp.-Start_Position.bp.)/Length.snp. < 1000000)

if (nrow(segments) > 0) 
	{
	for (i in 1:nrow(segments))
	{
		temp=NULL
		temp=as.data.frame(segments[i,c("Chromosome","Start_Position.bp.","End_Position.bp.")])
		colnames(temp)=c("chr","start","end")
		results=AnnotateSegment(temp)
		if (nrow(results) > 0) 
			{
			results=results[,c("chr", "start", "end", "geneName","geneID")]
			colnames(results)=c("Chrom", "Start", "End", "Gene","GeneID")
			results[,"Chrom"]=gsub("chr","",results[,"Chrom"])
			results$TITAN=as.data.frame(segments)[i,"TITAN_call"]
			results$Name=name
			results=results[,c("Name", "Chrom", "Start", "End", "TITAN","Gene","GeneID")]
			loh=rbind(loh,results)
			}
	}

	if (species=="Mouse")
	{
		segments.gr=makeGRangesFromDataFrame(segments,seqnames.field="Chromosome", start.field="Start_Position.bp.", end.field="End_Position.bp.", keep.extra.columns=T)
		ncruc.gr=GRanges(4, IRanges(89311040, 89511040))
		olaps=suppressWarnings(findOverlaps(segments.gr,ncruc.gr))
		ncruc=data.frame(pintersect(segments.gr[queryHits(olaps)], ncruc.gr[subjectHits(olaps)]))[,c("seqnames","start","end","TITAN_call")]
		if (nrow(ncruc) > 0)
		{
			ncruc$Gene="Cdkn2_ncruc"
			ncruc$Name=name
			ncruc$GeneID=NA
			colnames(ncruc)=c("Chrom", "Start", "End", "TITAN", "Gene", "Name","GeneID")
			ncruc=ncruc[,c("Name","Chrom", "Start", "End", "TITAN", "Gene","GeneID")]
			loh=rbind(loh,ncruc)
		}
	}
	loh=tbl_df(loh) %>% mutate_each(as.character)

	loh = loh %>% 
	#group_by(Gene) %>% 
	#arrange(TITAN,.by_group=T) %>% 
	#filter(row_number()==1) %>% 
	arrange(as.numeric(Chrom),as.numeric(Start),as.numeric(End))
} else {
		loh[1,"Name"]=name
		loh[1,"Chrom"]=1
		loh[1,"Start"]=1
		loh[1,"End"]=1
		loh[1,"TITAN"]="NLOH"
		loh[1,"Gene"]="EMPTY"
		loh[1,"GeneID"]="EMPTY"
}

write.table(loh,paste(name,"/results/LOH/",name,".LOH.genes.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

CGC=read.delim(CGC,header=T,sep="\t")

loh_cgc = loh %>%
filter(Gene %in% c(as.character(CGC[,1]),"Cdkn2_ncruc"))

write.table(loh_cgc,paste(name,"/results/LOH/",name,".LOH.genes.CGC.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")