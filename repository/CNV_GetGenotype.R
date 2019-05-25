#!/usr/bin/Rscript

##########################################################################################
##
## CNV_GetGenotype.R
##
## Extract genotypes for different mouse allels.
##
##########################################################################################

args <- commandArgs(TRUE)

name = args[1]
genecode = args[2]
transcript = args[3]
allele = args[4]
position = args[5]
wt_allele = args[6]
del_allele = args[7]

suppressMessages(library(ggplot2))

wt_allele=as.numeric(unlist(strsplit(wt_allele,",")))
del_allele=as.numeric(unlist(strsplit(del_allele,",")))

normal=read.delim(paste(name,"/results/Genotype/",name,".Normal.Genotypes.CNV.txt",sep=""))
tumor=read.delim(paste(name,"/results/Genotype/",name,".Tumor.Genotypes.CNV.txt",sep=""))

colnames(normal)=c("Chrom","Pos","Depth")
colnames(tumor)=c("Chrom","Pos","Depth")

all_exons=readRDS(genecode)

exons=all_exons[all_exons[,"transcriptID"]==transcript,c("exonNum","chr","start","end")]

exons[,"chr"]=gsub("chr","",exons[,"chr"])

colnames(exons)=c("Exon", "Chrom", "Start", "End")

results=data.frame(matrix(ncol = 3,nrow=nrow(exons)))

colnames(results)=c("Exon","Type","Median_depth")

results[,"Exon"]=as.character(results[,"Exon"])

for (i in 1:nrow(exons))
{
	results[i,"Type"]="Normal"
	results[i,"Exon"]=as.character(exons[i,"Exon"])
	results[i,"Median_depth"]=median(normal[normal[,"Pos"] > exons[i,"Start"] & normal[,"Pos"] < exons[i,"End"],"Depth"])
}

for (i in 1:nrow(exons))
{
	results[i+nrow(exons),"Type"]="Tumor"
	results[i+nrow(exons),"Exon"]=as.character(exons[i,"Exon"])
	results[i+nrow(exons),"Median_depth"]=median(tumor[tumor[,"Pos"] > exons[i,"Start"] & tumor[,"Pos"] < exons[i,"End"],"Depth"])
}

results[,"Exon"]=as.numeric(results[,"Exon"])

output=data.frame(matrix(ncol = 9,nrow=2))

colnames(output)=c("Name","Allele","CHROM","POS","REF","ALT","Count_Ref","Count_Alt","Comment")

output[1,"Count_Alt"]=median(results[del_allele,"Median_depth"])
output[1,"Count_Ref"]=median(results[wt_allele,"Median_depth"])
output[2,"Count_Alt"]=median(results[nrow(exons)+del_allele,"Median_depth"])
output[2,"Count_Ref"]=median(results[nrow(exons)+wt_allele,"Median_depth"])
output[1,"Name"]=paste(name,"_Normal",sep="")
output[2,"Name"]=paste(name,"_Tumor",sep="")
output[,"Allele"]=allele
output[,"CHROM"]=unlist(strsplit(position,":"))[1]
output[,"POS"]=unlist(strsplit(position,":"))[2]
output[,"REF"]=paste("Exon ",paste(wt_allele,sep="",collapse=","),sep="")
output[,"ALT"]=paste("Exon ",paste(del_allele,sep="",collapse=","),sep="")
output[,"Comment"]=""
write.table(output,paste(name,"/results/Genotype/",name,".Genotypes.temp.CNV.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

suppressWarnings(ggplot(data=results, aes(x=Exon,y=Median_depth,fill=factor(Type))) +
geom_bar(stat="identity",position="dodge") + theme(axis.line = element_line(colour = "black")) +
scale_fill_discrete(name="Type", labels=c("Tumor", "Normal")) +
xlab("Exons") + ylab("Median Depth") + theme_bw() + coord_cartesian(ylim=c(0,200)) +
ggtitle(paste("Mouse:    ",name,"Allele: ",allele,sep="\t")) + theme(plot.title = element_text(hjust = 0.5)) +
scale_x_continuous(breaks=1:ceiling(round((nrow(results)/2),digits=1))) +
scale_fill_manual("legend", values = c("Tumor" = "black", "Normal" = "grey")))

ggsave(paste(name,"/results/Genotype/",name,".",allele,".CNV.pdf",sep=""),width = 8, height = 4)

invisible(dev.off())