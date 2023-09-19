#!/usr/bin/env Rscript

# Script to automate launch script and input.tsv table creation for remapping
# see https://github.com/roland-rad-lab/MoCaSeq/tree/human-pipeline-nextflow-2
# @author: marcus.wagner@tum.de
# 
# This script generates runner scripts and input tables for remapping MoCaSeq nextflow runs on LRZ.
# Please make sure all directories you want to use are valid and you have permissions.
# 
# Cluster information:
# This software is intended to be executed on the CoolMUC-3 Linux cluster of LRZ.
# With 48h runtime samples will be processed in smaller batches.
# Assumed runtime: 0.1h per 1GB of data

# example content of input.tsv 
# Sample_Name	Sample_Group	Library_ID	Lane	Colour_Chemistry	SeqType	Organism	Type	R1	R2	BAM
# PCSI_0357_St_R	PCSI_0357	?	?	?	wgs	Human	Normal	NA	NA	/gpfs/scratch/pn29ya/ga89tog2/ga89tog2/test/results/GRCh38.p12/PCSI_0357_St_R/results/bam/PCSI_0357_St_R.Normal.bam
# PCSI_0357_Pa_P_5262	PCSI_0357	?	?	?	wgs	Human	Tumor	NA	NA	/gpfs/scratch/pn29ya/ga89tog2/ga89tog2/test/results/GRCh38.p12/PCSI_0357_Pa_P_5262/results/bam/PCSI_0357_Pa_P_5262.Tumor.bam

library(data.table)
library(dplyr)
library(tidyr)
library(splitstackshape)
library(readxl)

# MoCaSeq call parameters
# make sure these directories exist or adapt these values to your need
ccc_dir <- '${HOME}/images-live/mocaseq2'
project_dir <- '/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS'
work_dir <- '/gpfs/scratch/pn29ya/$USER/mocaseq-nextflow/remap'
reference_dir <- '/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/references/bashMoCaSeq'
repo_dir <- '/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/software/MoCaSeq'
genome_build.human <- 'GRCh38.p12'
custom_config_version <- 'serial-std'
custom_config_base <- file.path(project_dir, 'MoCaSeq', 'conf')
out_base <- file.path(project_dir, 'input')

# batch01 test samples for benchmarking to previous work
# batch01_test <- c("RAMP_0007_Lv_M_5262", "RAMP_0007_Mu_R", "RAMP_0008_Ln_M_526",
#                   "RAMP_0008_Lv_M_526", "RAMP_0008_Mu_R", "RAMP_0008_Pa_P")

# read master table
dt.master <- read_xlsx('data/data/EGA_mastertable.xlsx') %>% as.data.table()
dt.master[, SampleID := gsub('.bam', '', FileName)]

# get only samples that are currently present on LRZ
# alternatively filter for samples here
# TODO dt.master[BatchSub %in% c('4d', '4e', '4f') & !is.na(hg19bam)]
dt.remap <- dt.master[!is.na(hg19bam) & !remaped,
                      .(SampleID, Batch, DonorName, hg19bam, userREMAP)]

# input table columns
# Sample_Name, Sample_Group, Library_ID, Lane, Colour_Chemistry, SeqType, Organism, Type, R1, R2, BAM
dt.remap[, ':=' (Sample_Name = SampleID, Sample_Group = DonorName,
                 Library_ID = '?', Lane = '?', 
                 Colour_Chemistry = '?', SeqType = 'wgs', Organism = 'Human', 
                 BAM = NA)]

# set sample type by grabbing Standard/Reference identifier 'St_R'
dt.remap[grep('_R', SampleID), Type := 'Normal']
dt.remap[is.na(Type), Type := 'Tumor']

# create R1 and R2 reference on bam files
dt.remap[, ':=' (R1 = hg19bam, R2 = hg19bam, hg19bam = NULL)]

# remove unneccesary information
dt.remap[, ':=' (DonorName = NULL, SampleID = NULL)]

# get sample groups
sample_groups <- dt.remap[, paste0(Batch, '-', Sample_Group)] %>% unique()

# this function creates all variables specific for a certain sample run
sample_call <- function(dt){
paste0('
# specify sample
sample=',dt$Sample_Name,'
bamDir=',dirname(dt$R1),'
bamName=',basename(dt$R1),'
bamType="',dt$Type,'"
# submit subjob for sample remapping
srun --ntasks=1 --exclusive --mem 45568mb -J $sample -o ./', dt$Sample_Name,
'.%j.%N.out ${mocaseqDir}/launch/ccc_remap_wrapper.sh -ccc $cccDir -wd $workingDir -m $mocaseqDir -bd $bamDir -bf $bamName -s $sample -rd $referencesDir -t $bamType > ${sample}-remap.out & 
sleep 4')
}

# create these files for each sample group:
#  - runner script for bash MoCaSeq
#  - input table for .nextflow execution
#  - nextflow execution call as bash script
for (s_group in sample_groups) {
  # create launch dir for sample
  sample_remap_dir <- file.path(getwd(), 'launch', 'remap', s_group)
  system(paste0('mkdir -p ', sample_remap_dir))
  
#   # write caller script file by Sample_Group for bash MoCaSeq
#   cat(paste0('#!/bin/bash
# #SBATCH -J ', s_group, '-remap
# #SBATCH -D ./
# #SBATCH --clusters=mpp3
# #SBATCH --partition=mpp3_batch
# #SBATCH --nodes=1
# #SBATCH --time=48:00:00
# #SBATCH --ntasks-per-node=2
# #SBATCH --cpus-per-task=32
# #SBATCH --mem=91136mb
# #SBATCH --get-user-env
# #SBATCH --mail-type=end
# #SBATCH --mail-user=marcus.wagner@tum.de
# #SBATCH --export=NONE
# 
# # SLURM script for remapping sample group ', s_group, ' .bam files hg19 -> hg38.
# # Execution on a single node of CoolMUC-3 cluster.
# # ClusterInfo: mpp3_batch: max 48h, 90GB, 64 CPUs, jobs(50,dynamic)
# # See https://doku.lrz.de/job-processing-on-the-linux-cluster-10745970.html
# # !!! Please check:
# # - you have the mocaseq2 charldiecloud container in $HOME/images-live
# # - all paths are accessible, expecially the reference path 
# 
# module load slurm_setup
# 
# # set container path
# cccDir=${HOME}/images-live/mocaseq2
# 
# # specify mount paths
# workingDir=/gpfs/scratch/pn29ya/$USER/mocaseq-slurm/remap
# mocaseqDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/software/MoCaSeq/
# referencesDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/references/bashMoCaSeq/
# # make sure the working dir exists and create a symbolic link to the references named "ref" inside working directory
# mkdir -p $workingDir
# mkdir -p /gpfs/scratch/pn29ya/$USER/tmp
# ln -s $referencesDir ref
# mv ref $workingDir
# ', paste0(sample_call(dt.remap[Sample_Group == strsplit(s_group, '-')[[1]][2]]), collapse = '\n'),'
# 
# wait # for completion of background tasks
# '),
#       file = file.path(sample_remap_dir, paste0(s_group, '_cm3.sh')))
#   
#   # make runner file executable for user
#   system(paste0("chmod ug+x ",
#                 file.path(sample_remap_dir, paste0(s_group, '_cm3.sh'))))
  
  # .tsv file for nextflow remap
  write.table(dt.remap[Sample_Group == strsplit(s_group, '-')[[1]][2], 
                         .(Sample_Name, Sample_Group, Library_ID, Lane,
                           Colour_Chemistry, SeqType, Organism, Type, R1, R2, BAM)],
              file = file.path(sample_remap_dir, paste0(s_group, '.tsv')),
              sep = '\t', row.names = F, quote = F)
  
  # nextflow call
  # write caller script file by Sample_Group
  cat(paste0('# MoCaSeq launch script for sample group ', s_group, '\n'),
      paste0('projectDir=', project_dir, '\n'),
      paste0('workDir=', work_dir, '\n'),
      paste0('nextflow run ', file.path(repo_dir, 'main.nf'),
             ' -profile charliecloud,slurm', ' -entry MAP', ' -with-report -with-timeline',
             ' -work-dir ', file.path('${workDir}', 'remap', 'work'),
             ' --output_base ', file.path('${projectDir}', 'input'),
             ' --genome_build.human ', genome_build.human,
             ' --custom_config_version ', custom_config_version,
             ' --custom_config_base ', file.path(repo_dir, 'conf'),
             ' --input ', paste0(s_group, '.tsv\n')),
      file = file.path(sample_remap_dir, paste0(s_group, '_nf.sh')))
  
  # make runner file executable for user
  system(paste0("chmod ug+x ", 
                file.path(sample_remap_dir, paste0(s_group, '_nf.sh'))))
}