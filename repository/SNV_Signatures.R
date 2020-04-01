#!/usr/bin/Rscript

##########################################################################################
##
## SNV_Signatures.R
##
## Get canonical signatures for samples
##
##########################################################################################

args <- commandArgs(TRUE)

name = args[1]

library(SomaticSignatures)
library(SomaticCancerAlterations)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(datasets)
library(dplyr)
library(tidyr)
library(deconstructSigs)

file=paste0(name,"/results/Mutect2/",name,".Mutect2.vcf")
sampledf=data.frame()
t=read.table(file, header=F, sep="\t")
t=t[,c(1,2,4,5)]
colnames(t)=c("chr","pos","ref","alt")
t$Sample=name
t=t[,c("Sample","chr","pos","ref","alt")]
sampledf=t
sampledf=sampledf %>% filter(chr %in% c(1:19,"X","Y"))
sampledf$chr=paste0("chr",sampledf$chr)

sigs.input <- mut.to.sigs.input(mut.ref = sampledf, sample.id = "Sample", chr = "chr", pos = "pos", ref = "ref", alt = "alt", bsg=BSgenome.Mmusculus.UCSC.mm10)

sample = whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.nature2013, associated=c("Signature.1A", "Signature.2",  "Signature.3",  "Signature.4",  "Signature.5",  "Signature.6",  "Signature.7",  "Signature.8",  "Signature.9",  "Signature.10", "Signature.11", "Signature.12", "Signature.13", "Signature.14", "Signature.15", "Signature.16", "Signature.17", "Signature.18", "Signature.19", "Signature.20", "Signature.21"), sample.id = name, contexts.needed = TRUE, tri.counts.method = 'default', signature.cutoff = 0.2)
pdf(paste0(name,"/results/Mutect2/",name,"_Nature_Pie.pdf",sep=""))
makePie(sample)
dev.off()
pdf(paste0(name,"/results/Mutect2/",name,"_Nature_Bar.pdf",sep=""))
plotSignatures(sample)
dev.off()

sample = whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.cosmic, sample.id = name, contexts.needed = TRUE, tri.counts.method = 'default', signature.cutoff = 0.2)
pdf(paste0(name,"/results/Mutect2/",name,"_Cosmic_Pie.pdf",sep=""))
makePie(sample)
dev.off()
pdf(paste0(name,"/results/Mutect2/",name,"_Cosmic_Bar.pdf",sep=""))
plotSignatures(sample)
dev.off()