#!/usr/bin/Rscript

##########################################################################################
##
## SNV_Select",input_format,"Output.R
##
## Filters output for human annoted files.
##
##########################################################################################

args <- commandArgs(TRUE)

name = args[1]
input_format = args[2]
species= args[3]
CGC = args[4]
TruSight = args[5]

library(data.table)

file=as.data.frame(fread(paste(name,"/results/",input_format,"/",name,".",input_format,".txt",sep="")),header=T,sep="\t")

if (species=="Human")
{
	file$AF[is.na(file$AF)]=0
	file$AC[is.na(file$AC)]=0
	file$AN[is.na(file$AN)]=0

	sel=file[file[,"AF"] <0.1 & file[,"G5"]=="FALSE" & file[,"AN"] < 100 | file[,"AF"] <0.01 & file[,"AN"] >= 100 & file[,"G5"]=="FALSE",] 

	sel[is.na(sel)]=" "

	write.table(sel,paste(name,"/results/",input_format,"/",name,".",input_format,".NoCommonSNPs.txt",sep=""), col.names=T,row.names=F, quote=F, sep="\t")

	sel=sel[sel[,"ANN[*].IMPACT"] %in% c("HIGH", "MODERATE"),]

	sel[is.na(sel)]=" "

	write.table(sel,paste(name,"/results/",input_format,"/",name,".",input_format,".NoCommonSNPs.OnlyImpact.txt",sep=""), col.names=T,row.names=F, quote=F, sep="\t")

	CGC=read.delim(CGC,header=T,sep="\t")

	sel=sel[sel[,"ANN[*].GENE"]%in% CGC[,1],]

	sel[is.na(sel)]=" "

	write.table(sel,paste(name,"/results/",input_format,"/",name,".",input_format,".NoCommonSNPs.OnlyImpact.CGC.txt",sep=""), col.names=T,row.names=F, quote=F, sep="\t")

	TruSight=read.delim(TruSight,header=T,sep="\t")

	sel=sel[sel[,"ANN[*].GENE"]%in% TruSight[,1],]

	sel[is.na(sel)]=" "

	write.table(sel,paste(name,"/results/",input_format,"/",name,".",input_format,".NoCommonSNPs.OnlyImpact.TruSight.txt",sep=""), col.names=T,row.names=F, quote=F, sep="\t")
} else if (species=="Mouse") {
	sel=file[file[,"ANN[*].IMPACT"] %in% c("HIGH", "MODERATE"),]

	sel[is.na(sel)]=" "

	write.table(sel,paste(name,"/results/",input_format,"/",name,".",input_format,".NoCommonSNPs.OnlyImpact.txt",sep=""), col.names=T,row.names=F, quote=F, sep="\t")

	CGC=read.delim(CGC,header=T,sep="\t")

	sel=sel[sel[,"ANN[*].GENE"]%in% CGC[,1],]

	sel[is.na(sel)]=" "

	write.table(sel,paste(name,"/results/",input_format,"/",name,".",input_format,".NoCommonSNPs.OnlyImpact.CGC.txt",sep=""), col.names=T,row.names=F, quote=F, sep="\t")
}

