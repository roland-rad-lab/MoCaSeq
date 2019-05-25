#!/usr/bin/Rscript

##########################################################################################
##
## CNV_EstimateCoverage.R
##
## Estimate coverage, primarily used for lcWGS.
##
##########################################################################################

args <- commandArgs(TRUE)

name = args[1]

coverage=data.frame()

normal_metrics <- read.csv(paste0(name,"/results/QC/",name,".Normal.bam.metrics"),skip=10,header=T,sep="\t",stringsAsFactors=F)
tumor_metrics <- read.csv(paste0(name,"/results/QC/",name,".Tumor.bam.metrics"),skip=10,header=T,sep="\t",stringsAsFactors=F)

normal_mean <- sum(normal_metrics$coverage_or_base_quality*(normal_metrics$high_quality_coverage_count/sum(normal_metrics$high_quality_coverage_count)))
tumor_mean <- sum(tumor_metrics$coverage_or_base_quality*(tumor_metrics$high_quality_coverage_count/sum(tumor_metrics$high_quality_coverage_count)))
normal_dist = cumsum(normal_metrics[,2])/sum(normal_metrics[,2])
normal_median = which(abs(normal_dist-0.5)==min(abs(normal_dist-0.5)))
tumor_dist = cumsum(tumor_metrics[,2])/sum(tumor_metrics[,2])
tumor_median = which(abs(normal_dist-0.5)==min(abs(normal_dist-0.5)))

coverage = rbind(coverage,data.frame(name=name,Type="Normal",MeanCov=round(normal_mean,3),MedianCov=normal_median))
coverage = rbind(coverage,data.frame(name=name,Type="Tumor",MeanCov=round(tumor_mean,3),MedianCov=tumor_median))

write.table(coverage,paste0(name,"/results/QC/",name,".Normal.bam.coverage"),quote=F,row.names=F,col.names=T,sep="\t")

