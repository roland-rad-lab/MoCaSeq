#!/usr/bin/Rscript

##########################################################################################
##
## SV_SelectGenesFromManta.R
##
## Selects genes from Manta for use in oncoprints.
##
##########################################################################################

args <- commandArgs(TRUE)

name = args[1]

library(data.table)
library(dplyr)
library(tidyr)
library(readr)

tab=read_tsv(paste(name,"/results/Manta/",name,".Manta.txt",sep=""))

colnames(tab)=gsub("\\[\\*]","",colnames(tab))
colnames(tab)=gsub("\\[","",colnames(tab))
colnames(tab)=gsub("\\]","",colnames(tab))

tab = tab %>% 
mutate_if(is.integer, as.numeric) %>%
mutate(GENTumor.SR.AF=GENTumor.SR1/(GENTumor.SR1+GENTumor.SR0)) %>%
mutate(GENTumor.PR.AF=GENTumor.PR1/(GENTumor.PR1+GENTumor.PR0)) %>%
filter(GENTumor.SR1 + GENTumor.SR0 >= 10 | GENTumor.PR1 + GENTumor.PR0 >= 10) %>%
filter(ANN.IMPACT %in% c("HIGH","MODERATE")) %>%
filter(grepl("fusion|frameshift|stop|splice|start",ANN.EFFECT)) %>%
separate(ANN.GENE,into=c("Gene1","Gene2"),sep="&") %>% 
select(c("Gene1","Gene2","ANN.EFFECT")) %>%
gather(key,Gene,1:2) %>%
select(-key) %>%
filter(!is.na(Gene)) %>% 
separate(ANN.EFFECT,into=c("EFF1","EFF2"),sep="&") %>%
gather(key,Effect,1:2) %>% 
select(-key) %>%
filter(!is.na(Effect)) %>%
distinct()

write.table(tab,paste(name,"/results/Manta/",name,".Manta.genes.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")