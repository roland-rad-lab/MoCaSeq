#!/usr/bin/env Rscript

# Script to generate launch data (input table, nextflow and bash script) for MoCaSeq analysis
# see https://github.com/roland-rad-lab/MoCaSeq/tree/human-pipeline-nextflow-2
# @author: marcus.wagner@tum.de
# 
# This script generates launch scripts and input tables for MoCaSeq analyses executed on LRZ.
# Please make sure:
#  - to use the latest master table version.
#  - all directories you want to use are valid and you have permissions.
# 
# Example content of input table (use .tsv 
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
work_dir <- '/gpfs/scratch/pn29ya/$USER/mocaseq-nextflow/mocaseq'
reference_dir <- '/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/references/bashMoCaSeq'
repo_dir <- '/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/software/MoCaSeq'
genome_build.human <- 'GRCh38.p12'
custom_config_version <- 'serial-std'
custom_config_base <- file.path(project_dir, 'MoCaSeq', 'conf')
out_base <- file.path(project_dir, 'output')

# batch01 test samples for benchmarking to previous work
# batch01_test <- c("RAMP_0007_Lv_M_5262", "RAMP_0007_Mu_R", "RAMP_0008_Ln_M_526",
#                   "RAMP_0008_Lv_M_526", "RAMP_0008_Mu_R", "RAMP_0008_Pa_P")

# read master table
dt.master <- read_xlsx('data/data/EGA_mastertable.xlsx') %>% as.data.table()
dt.master[, SampleID := gsub('.bam', '', FileName)]

# get only samples that are currently present on LRZ
# alternatively filter for samples here
dt.mocaseq <- dt.master[ !is.na(hg38bam) & mocaseqed == F & remaped == T,
                        .(SampleID, Batch, DonorName, hg38bam, userMoCaSeq)]

# input table columns
# Sample_Name, Sample_Group, Library_ID, Lane, Colour_Chemistry, SeqType, Organism, Type, R1, R2, BAM
dt.mocaseq[, ':=' (Sample_Name = SampleID, Sample_Group = DonorName, Library_ID = '?', Lane = '?', 
                 Colour_Chemistry = '?', SeqType = 'wgs', Organism = 'Human', 
                 R1 = NA, R2 = NA, BAM = hg38bam,
                 hg38bam = NULL, DonorName = NULL)]


# set sample type by grabbing Standard/Reference identifier 'St_R'
dt.mocaseq[grep('_R', SampleID), Type := 'Normal']
dt.mocaseq[is.na(Type), Type := 'Tumor']

# remove unneccesary information
dt.mocaseq[, ':=' (SampleID = NULL)]

# get sample groups
sample_groups <- dt.mocaseq[, paste0(Batch, '-', Sample_Group)] %>% unique()

# --------- MoCaSeq Call loop ------------- 
# write input and runner files per sample group
for (s_group in sample_groups) {
  # create launch dir for sample !!! assumes to be executed in git repo root !!!
  sample_mocaseq_dir <- file.path(getwd(), 'launch', 'mocaseq', s_group)
  system(paste0('mkdir -p ', sample_mocaseq_dir))
  
  # write input.tsv files by Sample_Group
  write.table(dt.mocaseq[Sample_Group == strsplit(s_group, '-')[[1]][2], 
                         .(Sample_Name, Sample_Group, Library_ID, Lane,
                           Colour_Chemistry, SeqType, Organism, Type, R1, R2, BAM)],
            file = file.path(sample_mocaseq_dir, paste0(s_group, '.tsv')),
            sep = '\t', row.names = F, quote = F)
  
  # write caller script file by Sample_Group
  cat(paste0('# MoCaSeq launch script for sample group', s_group, '\n'),
      paste0('# path to DSS project folder\n','projectDir=', project_dir, '\n'),
      paste0('# path to MoCaSeq git repo\n', 'repoDir=', repo_dir, '\n'),
      paste0('# path to working dir\n', 'workDir=', work_dir, '\n'),
      paste0('# nextflow needs this variable to get SLURM job status\n',
             'export SLURM\\_CLUSTERS="', strsplit(custom_config_version, '-')[[1]][1], '"\n'),
      paste0('\n','nextflow run ', file.path('${repoDir}', 'main.nf'),
             ' -profile charliecloud,slurm -with-report -with-timeline',
             ' -work-dir ', '${workDir}', 
             ' --output_base ', file.path('${projectDir}', 'output'), 
             ' --genome_build.human ', genome_build.human,
             ' --custom_config_version ', custom_config_version,
             ' --custom_config_base ', file.path('${repoDir}', 'conf'), 
             ' -N marcus.wagner@tum.de --input ', paste0(s_group, '.tsv\n')),
      file = file.path(sample_mocaseq_dir, paste0(s_group, '_nf.sh')))
  
  # make runner file executable for user
  system(paste0("chmod ug+x ", 
                file.path(sample_mocaseq_dir, paste0(s_group, '_nf.sh'))))
}
