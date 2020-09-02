#!/usr/bin/Rscript

##########################################################################################
##
## CNV_RunFACETS.R
##
## Run FACETS.
##
##########################################################################################

#parallel 'Rscript '/media/rad/SSD1/DNA/repository/all_RunFACETS.R' {} Mouse 300 10' ::: mPDAC0004_7639_PPT-1

args <- commandArgs(TRUE)

name = args[1]
species = args[2]
cval <- args[3]
ndepth <- args[4]

suppressPackageStartupMessages(library(facets))
  
if (species=="Human"){
	gbuild="hg38"
} else if (species=="Mouse"){
	gbuild="mm10"
}

set.seed(1234)
datafile = paste0(name,"/results/FACETS/",name,".pileup.txt.gz")
rcmat = readSnpMatrix(datafile)
xx = preProcSample(rcmat, gbuild=gbuild,ndepth=ndepth)
oo=procSample(xx,cval=cval)
fit=suppressWarnings(emcncf(oo))

pdf(paste0(name,"/results/FACETS/",name,".FACETS.pdf"))
plotSample(x=oo,emfit=fit)
garbage <- dev.off()

pdf(paste0(name,"/results/FACETS/",name,".FACETS.spider.pdf"))
logRlogORspider(oo$out, oo$dipLogR)
garbage <- dev.off()

info=paste0("Purity: ",fit$purity," Ploidy: ",fit$ploidy)
write.table(info, paste0(name,"/results/FACETS/",name,".FACETS.txt"), col.names=F, row.names=F, quote=F, sep="\t")