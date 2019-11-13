# CALL THIS ONCE
library(data.table)
library(tidyr)
library(matrixStats)
library(doParallel)
registerDoParallel(cores=5)
getQC <- function(mousefolder){
  
  mouseid <- gsub("/", "", mousefolder)
  QCfolder <- paste0(mousefolder, "results/QC/", mouseid, "_data/")
  
  file1 <- paste0(QCfolder, "multiqc_general_stats.txt")
  file2 <- paste0(QCfolder, "mqc_picard_aligned_reads_1.txt")
  file3 <- paste0(QCfolder, "mqc_picard_percentage_target_bases_1.txt")
  
  file4 <- paste0(QCfolder, "mqc_picard_deduplication_1.txt") #mqc_picard_deduplication_1.txt (or multiqc_picard_dups.txt)
  if(!file.exists(file4)){
    file4 <- paste0(QCfolder, "multiqc_picard_dups.txt")
  }
  
  file5 <- paste0(QCfolder, "mqc_samtools-idxstats-xy-plot_1.txt")
  file6 <- paste0(QCfolder, "mqc_fastqc_per_sequence_quality_scores_plot_1.txt")
  file7 <- paste0(QCfolder, "mqc_fastqc_per_base_sequence_quality_plot_1.txt")
  file8 <- paste0(QCfolder, "multiqc_fastqc.txt")
  
  d1 <- fread(file1)
  d1 <- d1[, c("Sample","Picard_mqc-generalstats-picard-PERCENT_DUPLICATION")]
  names(d1) <- c("Sample", "Duplicates[%]")
  
  d2 <- fread(file2)
  names(d2) <- c("Sample", "AlignedReads", "UnalignedReads")
  
  d3 <- fread(file3, header=T)
  d3 <- d3[, c("Sample", "50")]
  names(d3) <- c("Sample", "MedianPercentageOfTargetBases(Coverage)")
  
  d4 <- fread(file4)
  d4[, rs := rowSums(d4[, -"Sample"])]
  d4data <- d4[, lapply(.SD, function(x){round(x / rs * 100, digits=2)}), .SDcols = !c("Sample", "rs")]
  d4 <- cbind(d4[, "Sample"], d4data)
  names(d4) <- gsub(" ", "",names(d4))
  names(d4) <- c("Sample", paste0(names(d4[, -c("Sample")]), "[%]"))
  
  d5 <- fread(file5)
  names(d5) <- c("Sample", "ChrY", "ChrX")
  
  d6 <- fread(file6, header=F, fill=T)
  for(row in 1:nrow(d6)){
    if(d6[row, V1] == ""){
      newname <- paste0(d6[row+1, V1], ".header")
      d6[row, V1 := newname]
    }
  }
  
  samples <- d6[!grepl(".header", V1), V1]
  
  buildD6 <- data.table()
  for(sample in samples){
    #sample <- "DS01_2259_LNMet-1.Normal.R1"
    
    valvec <- as.numeric(d6[V1 == sample, -c("V1")]) # get vector of count values
    maxvals <- sort(valvec, decreasing = T)[1:2] # find the 2 max values
    
    lastcol <- names(d6[, ncol(d6), with=F]) # find the name of the last column
    
    d6long <- data.table(gather(d6[V1 == sample], "column", "value", V2:lastcol)) # long to wide
    maxcolumns <- d6long[value %in% maxvals, column] # find colname with the max values
    
    phreadWithMaxScores <- d6[V1 == paste0(sample, ".header"), maxcolumns, with=F]
    names(phreadWithMaxScores) <- c("FastQCSequenceQuality[1]",	"FastQCSequenceQuality[2]")
    
    out <- data.table(Sample=sample, phreadWithMaxScores)
    buildD6 <- rbind(buildD6, out)
  }
  d6 <- copy(buildD6)
  
  d7 <- fread(file7, header=T, fill=T)
  # if colname Sample not found, rename V1
  if(!any(grepl("Sample",names(d7)))){
    setnames(d7, "V1", "Sample")
  }
  
  # also fix some broken files
  d7 <- d7[Sample != ""]
  
  rmean <- rowMeans(d7[, -c("Sample")])
  
  d7temp <- d7[, -c("Sample")]
  d7temp[, med := rowMedians(as.matrix(.SD))][]
  rmedian <- d7temp$med
  
  d7[, FastQCBaseSequenceQuality_Mean := rmean]
  d7[, FastQCBaseSequenceQuality_Median := rmedian]
  
  d7 <- d7[, .(Sample, FastQCBaseSequenceQuality_Mean, FastQCBaseSequenceQuality_Median)]
  
  d8 <- fread(file8)
  d8 <- d8[, .(Sample, overrepresented_sequences, adapter_content)]
  names(d8) <- c("Sample", "Overrepresented_sequences", "Adapter_content")
  
  mymerge <- function(x,y) merge(x,y,all=TRUE)
  out <- Reduce(mymerge,list(d1,d2,d3,d4,d5,d6,d7,d8))
  
  
  
  
  path9 <- paste0(mousefolder, "results/Copywriter/")
  files9 <- list.files(path9, full.names = T)
  files9 <- files9[grepl("pdf", files9)]
  files9 <- data.table(file.info(files9), keep.rownames = T) #bytes
  files9[, mb := size / 1000000]
  d9.minMB <- round(files9[, min(mb)], digits=2)
  
  file10 <- paste0(mousefolder, "results/Mutect2/",mouseid,".Mutect2.txt")
  d10.MB <- round(file.info(file10)$size / 1000000, digits=2)
  
  file11 <- paste0(mousefolder, "results/Mutect2/",mouseid,".m2.bam")
  d11.MB <- round(file.info(file11)$size / 1000000, digits=2)
  
  path12 <- paste0(mousefolder, "results/LOH/")
  files12 <- list.files(path12, full.names = T)
  files12 <- files12[grepl("pdf", files12)]
  files12 <- data.table(file.info(files12), keep.rownames = T) #bytes
  files12[, mb := size / 1000000]
  d12.minMB <- round(files12[, min(mb)], digits=2)
  
  file13 <- paste0(mousefolder, "results/LOH/",mouseid,".VariantsForLOH.txt")
  d13.MB <- round(file.info(file13)$size / 1000000, digits=2)
  
  out[, "CopywriteR[minMB]" := d9.minMB]
  out[, "Mutect2Calls[MB]" := d10.MB]
  out[, "Mutect2Bam[MB]" := d11.MB]
  out[, "LOHChrPDF[minMB]" := d12.minMB]
  out[, "LOHVariantsForLOH[MB]" := d13.MB]
  
  return(out)
}
#wdpath <- "/run/user/1000/gvfs/smb-share:server=imostorage.med.tum.de,share=fastq/Studies/AGRad_mPDAC/" #kata
wdpath <- "Y:/Studies/AGRad_mPDAC"
setwd(wdpath)


# get QC summary for 1 mouse
mousefolder <- "DS01_2259_LNMet-1/"
qc <- getQC(mousefolder)
qc # look at data in R
fwrite(output, "irgendwohin.txt", sep="\t") # write to file



# get QC summary for all mice in folder (~1 second per directory)
#dirs <- system("ls -d -- */", intern = T) # this works on linux
dirs <- list.dirs(recursive = F)
dirs <- dirs[!grepl("plot|annotated|misc|MultiSampleCalling|CohortDatabase", dirs)] # remove some other folders
dirs <- gsub("./", "", dirs)
dirs <- paste0(dirs, "/")

output <- foreach(i=1:length(dirs), .packages=c("data.table"), .combine=rbind) %dopar% {
  setwd(wdpath)
  library(tidyr)
  library(data.table)
  library(matrixStats)
  mousefolder <- dirs[i]
  qc <- getQC(mousefolder)
  qc
}
output #look at data
fwrite(output, "irgendwohin.txt", sep="\t") # write it to file
