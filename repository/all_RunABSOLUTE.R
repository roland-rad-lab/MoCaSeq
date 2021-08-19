#!/usr/bin/Rscript

##########################################################################################
##
## all_RunABSOLUTE.R
##
## Run ABSOLUTE.
##
##########################################################################################

#parallel 'Rscript /media/rad/SSD1/DNA/repository/all_RunABSOLUTE.R {} Copywriter Human' ::: $(ls)

args <- commandArgs(TRUE)

name = args[1]
cnv_method = args[2]
species = args[3]
mutect2 = args[4]
runmode = args[5] # SS or MS

suppressPackageStartupMessages(library(ABSOLUTE))

if (cnv_method=="Copywriter"){
	load(paste0(name,"/results/Copywriter/CNAprofiles/segment.Rdata"))
	segments = segment.CNA.object$output
	Selection = unique(grep("Normal",grep("Tumor",segments$ID,value=T),value=T))
	segments = segments[segments$ID==Selection,]
	segments = segments[segments$chrom<=22,]
	segments$loc.start = floor(segments$loc.start)
	segments$loc.end = floor(segments$loc.end)
	segments = segments[,c("chrom","loc.start","loc.end","num.mark","seg.mean")]
	colnames(segments)=c("Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
	write.table(segments, paste0(name,"/results/ABSOLUTE/",name,".Copywriter.seg.dat.fn"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	seg.dat.fn <- file.path(paste0(name,"/results/ABSOLUTE/",name,".Copywriter.seg.dat.fn"))
} else if (cnv_method=="HMMCopy") {
	segments = read.delim(paste0(name,"/results/HMMCopy/",name,".HMMCopy.20000.segments.txt"))
	segments = segments[segments[,"Chrom"] %in% seq(1,22,1),]
	if(species=="Human"){
	  chrom.sizes = c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,
	                  138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,
	                  83257441,80373285,58617616,64444167,46709983,50818468)
	  names(chrom.sizes) = c(1:22)
	}
	if(species=="Mouse"){
	  chrom.sizes = c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,
	                  130694993,122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566)
	  names(chrom.sizes) = c(1:19)
	}
	for (i in 1:nrow(segments)){
			segments[i,"Num_Probes"]=round((segments[i,"End"]-segments[i,"Start"])/sum(chrom.sizes)*1000000)
	}
	segments=segments[,c("Chrom","Start","End","Num_Probes","Mean")]
	colnames(segments)=c("Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
	write.table(segments, paste0(name,"/results/ABSOLUTE/",name,".HMMCopy.20000.segments.txt.fn"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	seg.dat.fn <- file.path(paste0(name,"/results/ABSOLUTE/",name,".HMMCopy.20000.segments.txt.fn"))
} else if(cnv_methood == "CNVKit"){
  
  if(runmode == "MS"){
    segfile <- paste0(name,"/results/CNVKit/matched/",name,".cns")
  } else if(runmode == "SS"){
    segfile <- paste0(name,"/results/CNVKit/single/",name,".Tumor.cns")
  } else {
    stop(paste0("Incorrect runmode parameter: ", runmode))
  }
  
  segments = read.delim(segfile)
  segments = segments[segments[,"chromosome"] %in% seq(1,22,1),]
  segments=segments[,c("chromosome","start","end","probes","log2")]
  colnames(segments)=c("Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
  write.table(segments, paste0(name,"/results/ABSOLUTE/",name,".CNVKit.seg.dat.fn"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
  seg.dat.fn <- file.path(paste0(name,"/results/ABSOLUTE/",name,".CNVKit.seg.dat.fn"))
}

# depending on if Mutect2 was called this is changed
if(mutect2 == "yes"){
  
  
  if(runmode == "MS"){
    maf=read.delim(paste0(name,"/results/Mutect2/",name,".Mutect2.vep.maf"),fill=T,skip=1)
    names(maf)[names(maf) == 't_maf'] <- 't_vaf'
    names(maf)[names(maf) == 'Start_Position'] <- 'Start_position'
    names(maf)[names(maf) == 'End_Position'] <- 'End_position'
    write.table(maf,paste0(name,"/results/Mutect2/",name,".Mutect2.vep.maf.fn"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    maf.fn <- file.path(paste0(name,"/results/Mutect2/",name,".Mutect2.vep.maf.fn"))
  } else if(runmode == "SS"){

    # identify if single sample is Tumor or Normal
    mafT=paste0(name,"/results/Mutect2/",name,".Tumor.Mutect2.vep.maf")
    mafN=paste0(name,"/results/Mutect2/",name,".Normal.Mutect2.vep.maf")
    
    if(file.exists(mafT)){
      ssPrefix <- "Tumor"
    } else if(file.exists(mafN)){
      ssPrefix <- "Normal"
    } else {
      stop("No MAF file found!") 
    }
    
    maf=read.delim(paste0(name,"/results/Mutect2/",name,".",ssPrefix,".Mutect2.vep.maf"),fill=T,skip=1)
    names(maf)[names(maf) == 't_maf'] <- 't_vaf'
    names(maf)[names(maf) == 'Start_Position'] <- 'Start_position'
    names(maf)[names(maf) == 'End_Position'] <- 'End_position'
    write.table(maf,paste0(name,"/results/Mutect2/",name,".",ssPrefix,".Mutect2.vep.maf.fn"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    maf.fn <- file.path(paste0(name,"/results/Mutect2/",name,".",ssPrefix,".Mutect2.vep.maf.fn"))
    
  } else {
    stop(paste0("Incorrect runmode parameter: ", runmode))
  }
  
  min.mut.af <- 0
} else if(mutect2 == "no"){
  maf.fn <- NULL
  min.mut.af <- NULL
} else {
  stop(paste0("Unknown value for argument 4 (",mutect2,")! This value of Mutect2 should be either yes or no"))
}

# DEFAULT PARAMETERS: https://www.genepattern.org/modules/docs/ABSOLUTE
sigma.p <- 0
max.sigma.h <- 0.015
min.ploidy <- 0.95
max.ploidy <- 10
max.as.seg.count <- 1500
max.non.clonal <- 0.05
max.neg.genome <- 0
copy_num_type <- "total"
primary.disease <- "NA"
platform <- "Illumina_WES"

if (species == "Mouse") genome_build="mm10" else if (species == "Human") genome_build="hg38"
genome_build <- genome_build
sample.name <- name
results.dir <- file.path(paste0(name,"/results/ABSOLUTE/"))

RunAbsolute(seg.dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy, primary.disease, platform, 
            sample.name, results.dir, max.as.seg.count, max.non.clonal, max.neg.genome, copy_num_type, 
            maf.fn=maf.fn, min.mut.af=min.mut.af, output.fn.base=NULL, verbose=T)


# read the results and print it like a normal human being
abs.rawfile <- paste0(results.dir, "/", name, ".ABSOLUTE.RData")
load(abs.rawfile)
seg.list <- list(seg.dat)
seg.results <- seg.list[[1]][["mode.res"]][["mode.tab"]]

# ABSOLUTE sometimes just can not process samples (error message: "Sample has failed ABSOLUTE")
abs.tabfile <- paste0(results.dir, "/", name, ".ABSOLUTE.results.tsv")
if(is.null(seg.results)){
  fileConn<-file(abs.tabfile)
  writeLines(c("Sample has failed ABSOLUTE"), fileConn)
  close(fileConn)
} else {
  colnames(seg.results)[1] <- "purity"
  colnames(seg.results)[2] <- "ploidy"
  write.table(seg.results, abs.tabfile, quote = F, sep="\t", row.names = F)
}


