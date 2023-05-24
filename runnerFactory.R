#!/usr/bin/env Rscript

# script to automate input.tsv table creation for mocaseq pipeline run
# @author: marcus.wagner@tum.de

# example content of input.tsv 
# Sample_Name	Sample_Group	Library_ID	Lane	Colour_Chemistry	SeqType	Organism	Type	R1	R2	BAM
# PCSI_0357_St_R	PCSI_0357	?	?	?	wgs	Human	Normal	NA	NA	/gpfs/scratch/pn29ya/ga89tog2/ga89tog2/test/results/GRCh38.p12/PCSI_0357_St_R/results/bam/PCSI_0357_St_R.Normal.bam
# PCSI_0357_Pa_P_5262	PCSI_0357	?	?	?	wgs	Human	Tumor	NA	NA	/gpfs/scratch/pn29ya/ga89tog2/ga89tog2/test/results/GRCh38.p12/PCSI_0357_Pa_P_5262/results/bam/PCSI_0357_Pa_P_5262.Tumor.bam

library(data.table)
library(dplyr)
library(tidyr)
library(splitstackshape)

dt.info <- fread('../../AllCompassAssociatedEgadFiles_2022-02-24.csv')
dt.info[Downloaded == 'corrupted-2021-12', Downloaded := 'corrupted_2021-12']
# separate batch
dt.info <- separate(dt.info, col = 'Downloaded', into = c('Batch', 'Date'), sep = '_') %>% 
  as.data.table()

# set column for EGA source data availability
dt.info[, EGAonLRZ := F]
dt.info[grep('yes', Date), EGAonLRZ := T]

# split yes/no from Date
dt.info[, Date := substr(Date, 1, 7)]

# MoCaSeq call parameters
# make sure these directories exist in your $USER paths
ega_data_path <- '/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/projects/hPDAC/ICGC_PACA_CA_WGS/input/GRCh37_bam'
project_dir <- '/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS'
work_dir <- '/gpfs/scratch/pn29ya/${USER}/${USER}/mocaseq-nf' # add [remap|mocaseq]/work later
genome_build.human <- 'GRCh38.p12'
custom_config_version <- 'mocaseq-lrz'
custom_config_base <- 'https://raw.githubusercontent.com/roland-rad-lab/MoCaSeq/human-pipeline-nextflow-2/conf'
input_prefix <- '/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/mocaseq_input_tables'

remap_out_base <- file.path(project_dir, 'input')
mocaseq_out_base <- file.path(project_dir, 'output')

# batch01 test samples
# batch01_test <- c("RAMP_0007_Lv_M_5262", "RAMP_0007_Mu_R", "RAMP_0008_Ln_M_526",
#                   "RAMP_0008_Lv_M_526", "RAMP_0008_Mu_R", "RAMP_0008_Pa_P")

# get only samples that are currently present on LRZ
# alternatively filter for samples here
dt.remap <- dt.info[EGAonLRZ == TRUE, .(FileID, FileName, SampleID, Batch)]

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
            file = file.path('input', 'remap', paste0(s_group, '.tsv')),
            sep = '\t', row.names = F, quote = F)

  # write caller script file by Sample_Group
  cat(paste0('nextflow run roland-rad-lab/MoCaSeq -r human-pipeline-nextflow-2 ',
             ' -profile charliecloud,slurm', ' -entry MAP',
             ' -work-dir ', file.path(work_dir, 'remap', 'work'),
             ' --output_base ', remap_out_base,
             ' --genome_build.human ', genome_build.human,
             ' --custom_config_version ', custom_config_version,
             ' --custom_config_base ', custom_config_base,
             ' --input ', input_prefix, 'remap/', s_group, '.tsv'),
      file = file.path('launch', 'remap', paste0(s_group, '_remap.sh')))
  
  # make runner file executable for user
  system(paste0("chmod u+x ",
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
            file = file.path('input', 'mocaseq', paste0(s_group, '.tsv')),
            sep = '\t', row.names = F, quote = F)
  
  # write caller script file by Sample_Group
  cat(paste0('nextflow run roland-rad-lab/MoCaSeq -r human-pipeline-nextflow-2 ',
             ' -profile charliecloud,slurm',
             ' -work-dir ', file.path(work_dir, 'mocaseq', 'work'), 
             ' --output_base ', mocaseq_out_base, 
             ' --genome_build.human ', genome_build.human,
             ' --custom_config_version ', custom_config_version,
             ' --custom_config_base ', custom_config_base, 
             ' --input ', input_prefix, 'mocaseq/', paste0(s_group, '.tsv')),
      file = file.path('launch', 'mocaseq', paste0(s_group, '_mocaseq.sh')))
  
  # make runner file executable for user
  system(paste0("chmod u+x ", 
                file.path('launch', 'mocaseq', paste0(s_group, '_mocaseq.sh'))))
}
