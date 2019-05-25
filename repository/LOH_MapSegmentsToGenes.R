#!/usr/bin/Rscript

##########################################################################################
##
## LOH_MapSegmentsToGenes.R
##
## Takes the segment file from either HMMCopy or Copywriter and maps them to genes.
##
##########################################################################################

#parallel 'Rscript /media/rad/SSD1/DNA/repository/LOH_MapSegmentsToGenes.R {} Mouse' ::: $(ls)

args <- commandArgs(TRUE)

name = args[1]
species = args[2]

suppressMessages(library(biomaRt))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(GenomicRanges))

if (species=="Mouse")
{
	dataset="mmusculus_gene_ensembl"
} else if (species=="Human") {
	dataset="hsapiens_gene_ensembl"
}

mart=useMart(biomart = 'ensembl', dataset = dataset, host = "useast.ensembl.org")

segments=paste(name,"/results/Titan/run_ploidy2/",name,"_cluster01.segs.txt",sep="")
loh=data.frame(Name=NULL,Chrom=NULL, Start=NULL, End=NULL, Mean=NULL,Gene=NULL)
segments=tbl_df(read.delim(segments))

segments = segments %>% filter(TITAN_call %in% c("ALOH","NLOH","DLOH","HOMD")) %>% filter(Length.snp. >= 10) %>% filter((End_Position.bp.-Start_Position.bp.)/Length.snp. < 1000000)

for (i in 1:nrow(segments))
{
	#print(paste("Getting information from segment ",i,"/",nrow(segments),".",sep=""))
	results=getBM(attributes=c("chromosome_name","start_position","end_position","external_gene_name"), 
		filters=c("chromosomal_region"),
		values=paste(segments[i,"Chromosome"],":",segments[i,"Start_Position.bp."],":",segments[i,"End_Position.bp."],sep=""),
		mart=mart)
	if (nrow(results) > 0) 
		{
		colnames(results)=c("Chrom", "Start", "End", "Gene")
		results$TITAN=as.data.frame(segments)[i,"TITAN_call"]
		results$Name=name
		results=results[,c("Name", "Chrom", "Start", "End", "TITAN","Gene")]
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
		colnames(ncruc)=c("Chrom", "Start", "End", "TITAN", "Gene", "Name")
		ncruc=ncruc[,c("Name","Chrom", "Start", "End", "TITAN", "Gene")]
		loh=rbind(loh,ncruc)
	}
}
loh=tbl_df(loh) %>% mutate_each(as.character)

loh = loh %>% 
#group_by(Gene) %>% 
#arrange(TITAN,.by_group=T) %>% 
#filter(row_number()==1) %>% 
arrange(as.numeric(Chrom),as.numeric(Start),as.numeric(End))

write.table(loh,paste(name,"/results/LOH/",name,".LOH.genes.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")