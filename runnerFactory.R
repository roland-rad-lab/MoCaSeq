#!/usr/bin/env Rscript

# Script to automate launch script and input.tsv table creation for MoCaSeq pipeline
# see https://github.com/roland-rad-lab/MoCaSeq/tree/human-pipeline-nextflow-2
# @author: marcus.wagner@tum.de
# 
# This script generates runner scripts and input tables for MoCaSeq runs on LRZ
# Please make sure all directories you want to use are valid and you have permissions.
# 

# example content of input.tsv 
# Sample_Name	Sample_Group	Library_ID	Lane	Colour_Chemistry	SeqType	Organism	Type	R1	R2	BAM
# PCSI_0357_St_R	PCSI_0357	?	?	?	wgs	Human	Normal	NA	NA	/gpfs/scratch/pn29ya/ga89tog2/ga89tog2/test/results/GRCh38.p12/PCSI_0357_St_R/results/bam/PCSI_0357_St_R.Normal.bam
# PCSI_0357_Pa_P_5262	PCSI_0357	?	?	?	wgs	Human	Tumor	NA	NA	/gpfs/scratch/pn29ya/ga89tog2/ga89tog2/test/results/GRCh38.p12/PCSI_0357_Pa_P_5262/results/bam/PCSI_0357_Pa_P_5262.Tumor.bam

library(data.table)
library(dplyr)
library(tidyr)
library(splitstackshape)

# MoCaSeq call parameters
# make sure these directories exist or adapt these values to your need
project_dir <- '/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS'
work_dir <- '/gpfs/scratch/pn29ya/ga89tog2/mocaseq-nextflow' # add [remap|mocaseq]/work later

ega_data_path <- '/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/projects/hPDAC/ICGC_PACA_CA_WGS/input/GRCh37_bam'
genome_build.human <- 'GRCh38.p12'
custom_config_version <- 'mocaseq-lrz'
custom_config_base <- file.path(project_dir, 'MoCaSeq', 'conf')
input_prefix <- file.path(project_dir, 'MoCaSeq', 'input')

remap_out_base <- file.path(project_dir, 'input')
mocaseq_out_base <- file.path(project_dir, 'output')

# batch01 test samples for benchmarking to previous work
# batch01_test <- c("RAMP_0007_Lv_M_5262", "RAMP_0007_Mu_R", "RAMP_0008_Ln_M_526",
#                   "RAMP_0008_Lv_M_526", "RAMP_0008_Mu_R", "RAMP_0008_Pa_P")

# read sample info table
dt.info <- fread('../../AllCompassAssociatedEgadFiles_2022-02-24.csv')
dt.info[Downloaded == 'corrupted-2021-12', Downloaded := 'corrupted_2021-12']
# separate batch
dt.info <- separate(dt.info, col = 'Downloaded', into = c('Batch', 'Date'), sep = '_') %>% 
  as.data.table()

# make a column for the sample status
# Possible StatusLRZ values: 'Not downloaded', 'EGA (hg19) bams available', 'Remap running',  'hg38 bam available', 'MoCaSeq running', 'MoCaSeq done'
dt.info[, StatusLRZ := 'Not downloaded']
dt.info[grep('yes', Date), StatusLRZ := 'EGA (hg19) bams available']

# split yes/no from Date
dt.info[, Date := substr(Date, 1, 7)]

# get only samples that are currently present on LRZ
# alternatively filter for samples here
dt.remap <- dt.info[StatusLRZ == 'EGA (hg19) bams available', .(FileID, FileName, SampleID, Batch)]

# input table columns
# Sample_Name, Sample_Group, Library_ID, Lane, Colour_Chemistry, SeqType, Organism, Type, R1, R2, BAM
dt.remap[, ':=' (Sample_Name = SampleID, Library_ID = '?', Lane = '?', 
                 Colour_Chemistry = '?', SeqType = 'wgs', Organism = 'Human', 
                 BAM = NA)]

# set sample type by grabbing Standard/Reference identifier 'St_R'
dt.remap[grep('_R', SampleID), Type := 'Normal']
dt.remap[is.na(Type), Type := 'Tumor']

# create R1 and R2 reference on bam files
dt.remap[, R1 := file.path(ega_data_path, Batch, FileID, FileName)]
dt.remap[, R2 := R1]

# extract Sample_Group information from SampleID
dt.remap <- separate(dt.remap, 'SampleID', into = c('Sample_Group', 'Tag'), sep = '_[A-Z]') %>%
  as.data.table()

# remove unneccesary information
dt.remap[, ':=' (FileID = NULL, FileName = NULL, Tag = NULL)]

# get sample groups
sample_groups <- dt.remap[, paste0(Batch, '-', Sample_Group)] %>% unique()

# write input and runner files per sample group for remapping
for (s_group in sample_groups) {
  # write input.tsv files by Sample_Group and with correct order of columns
  write.table(dt.remap[Sample_Group == strsplit(s_group, '-')[[1]][2], 
                       .(Sample_Name, Sample_Group, Library_ID, Lane,
                         Colour_Chemistry, SeqType, Organism, Type, R1, R2, BAM)],
            file = file.path('launch', 'remap', paste0(s_group, '.tsv')),
            sep = '\t', row.names = F, quote = F)
  
  # write caller script file by Sample_Group
  cat(paste0('# MoCaSeq remapping launch script for sample group', s_group, '\n'),
      paste0('projectDir=', project_dir, '\n'),
      paste0('workDir=', work_dir, '\n'),
      paste0('nextflow run ', file.path('${projectDir}', 'MoCaSeq', 'main.nf'),
             ' -profile charliecloud,slurm', ' -entry MAP',
             ' -work-dir ', file.path('${workDir}', 'remap', 'work'),
             ' --output_base ', file.path('${projectDir}', 'input'),
             ' --genome_build.human ', genome_build.human,
             ' --custom_config_version ', custom_config_version,
             ' --custom_config_base ', file.path('${projectDir}', 'MoCaSeq', 'conf'),
             ' --input ', file.path('${projectDir}', 'MoCaSeq', 'input', 'remap', s_group, '.tsv')),
      file = file.path('launch', 'remap', paste0(s_group, '_remap.sh')))
  
  # make runner file executable for user
  system(paste0("chmod ug+x ",
                file.path('launch', 'remap', paste0(s_group, '_remap.sh'))))
}

# --------- MoCaSeq Call ------------- 

# create input table for mocaseq call
dt.mocaseq <- copy(dt.remap)
dt.mocaseq[, ':=' (R2 = NA, R1 = NA,
                   BAM = file.path(remap_out_base, paste0(genome_build.human, '_bam'), 
                                   Batch, Sample_Name, 'results', 'bam', 
                                   paste0(Sample_Name, '.', Type, '.bam')))]

# write input and runner files per sample group
for (s_group in sample_groups) {
  # write input.tsv files by Sample_Group
  write.table(dt.mocaseq[Sample_Group == strsplit(s_group, '-')[[1]][2], 
                         .(Sample_Name, Sample_Group, Library_ID, Lane,
                           Colour_Chemistry, SeqType, Organism, Type, R1, R2, BAM)],
            file = file.path('launch', 'mocaseq', paste0(s_group, '.tsv')),
            sep = '\t', row.names = F, quote = F)
  
  # write caller script file by Sample_Group
  cat(paste0('# MoCaSeq launch script for sample group', s_group, '\n'),
      paste0('projectDir=', project_dir, '\n'),
      paste0('workDir=', work_dir, '\n'),
      paste0('nextflow run ', file.path('${projectDir}', 'MoCaSeq', 'main.nf'),
             ' -profile charliecloud,slurm',
             ' -work-dir ', file.path('${workDir}', 'mocaseq', 'work'), 
             ' --output_base ', file.path('${projectDir}', 'output'), 
             ' --genome_build.human ', genome_build.human,
             ' --custom_config_version ', custom_config_version,
             ' --custom_config_base ', file.path('${projectDir}', 'MoCaSeq', 'conf'), 
             ' --input ', file.path('${projectDir}', 'MoCaSeq', 'input', 'mocaseq', s_group, '.tsv')),
      file = file.path('launch', 'mocaseq', paste0(s_group, '_mocaseq.sh')))
  
  # make runner file executable for user
  system(paste0("chmod ug+x ", 
                file.path('launch', 'mocaseq', paste0(s_group, '_mocaseq.sh'))))
}
