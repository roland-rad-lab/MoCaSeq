# aspera content merge with AllEgaFiles table to master table

library(data.table)
library(readxl)
library(dplyr)
library(tidyr)
library(splitstackshape)

# read master table from Niklas
dt.master <- fread('EGA_mastertable.txt')

# read AllEgaFiles table
dt.info <- fread('../../AllCompassAssociatedEgadFiles_2022-02-24.csv')
dt.info[Downloaded == 'corrupted-2021-12', Downloaded := 'corrupted_2021-12']
# separate batch
dt.info <- separate(dt.info, col = 'Downloaded', into = c('Batch', 'Date'), sep = '_') %>% 
  as.data.table()

# set column for EGA source data availability
# Possible StatusLRZ values: 'Not downloaded', 'EGA (hg19) bams available', 'Remap running',  'hg38 bam available', 'MoCaSeq running', 'MoCaSeq done'
dt.info[, StatusLRZ := 'Not downloaded']
dt.info[grep('yes', Date), StatusLRZ := 'EGA (hg19) bams available']

# split yes/no from Date
dt.info[, Date := substr(Date, 1, 7)]

# investigare double SampleIDs
keep_samples <- dt.info[, .N, by = SampleID][N > 1, SampleID]
# library(ggplot2)
# ggplot(dt.info[SampleID %in% keep_samples[30:59]], aes(x = FileID, y = as.numeric(Bytes))) + geom_bar(stat = 'identity') + facet_wrap(~SampleID, scales = 'free')

# with double SampleIDs filter small files using IgnoreFile == T
dt.info[, IgnoreFile := F]
dt.info[SampleID %in% keep_samples & Bytes < 5e+10, IgnoreFile := T]

# read aspera content and split file path
dt.aspera <- fread('aspera_content.txt')
dt.aspera <- separate(dt.aspera, col = 'AsperaPath', into = c('EGAD', 'Part', 'Filename'), sep = '/') %>% as.data.table()
# dt.aspera <- separate(dt.aspera, col = 'Filename', sep = '[.]', into = c('EGAZ_SampleID', 'bam', 'c4gh')) %>% as.data.table()

# separate EGAZ id, data tags and Sample id
# EGAZ id
dt.aspera[, EGAZ := substr(Filename, 1, 15)]
# data tags
dt.aspera[, Tags := Filename]
dt.aspera[, Tags := gsub('EGAZ[0-9]{11}_', '', Tags)]
dt.aspera[, Tags := gsub('[A-Z]{4}[A-Z]?_[0-9]{4}_[A-Z][a-z]_[A-Z]_?([0-9]+)?([a-z]+)?', '', Tags)]
dt.aspera[, Tags := gsub('.bam.c4gh', '', Tags)]
# sample id
dt.aspera[, FileName := Filename]
dt.aspera[, FileName := gsub('EGAZ[0-9]{11}_', '', FileName)]
dt.aspera[, FileName := gsub('COMPASS_wgs_', '', FileName)]
dt.aspera[, FileName := gsub('COMPASS_wgs_COMPASS_Feasibility_1709_', '', FileName)]
dt.aspera[, FileName := gsub('PCSI_wgs_bam_head', '', FileName)]
dt.aspera[, FileName := gsub('PCSI_wgs_bam_', '', FileName)]
dt.aspera[, FileName := gsub('COMPASS_Feasibility_1709_', '', FileName)]
dt.aspera[, FileName := gsub('.c4gh', '', FileName)]

# remove unused columns
dt.aspera[, ':=' (Filename = NULL)] # , bam = NULL, c4gh = NULL

# sanity check to recreate file paths from dt.aspera
# dt.test <- dt.aspera[, file.path(EGAD, Part, paste0(EGAZ, '_', Tags, SampleID, '.bam.c4gh'))]
# write.table(dt.test, 'aspera_sanity_check.txt', quote = F, row.names = F, col.names = F)

# merge dt.aspera and dt.info into a master table
# dt.master <- merge(dt.info, dt.aspera, by = 'SampleID', all.x = T)

# 1. filter all small samples (Bytes < 5e+10) using IgnoreFile
# 2. filter all samples that were already processed (AnalysisFolderLRZ != 'yes-AllFiles')
# 3. filter all samples that are not observed in the aspera box (!is.na(EGAZ))
# 4. only download files that are not yet on LRZ StatusLRZ == 'Not downloaded'
# dt.filtered <- dt.master[IgnoreFile == F & StatusLRZ == 'Not downloaded' & AnalysisFolderLRZ == 'no' & !is.na(EGAZ)]
# setorder(dt.filtered, Batch)
# dt.filtered[, file.path(EGAD, Part, paste0(EGAZ, '_', Tags, SampleID, '.bam.c4gh'))] %>%
#   write.table(., file = 'aspera.filtered.sorted', quote = F, row.names = F, col.names = F)

# TODO: samples where EGAZ id is NA also have <SampleID>_(ref|tail|head) in aspera box
# EGAD00001004551/PART_2/EGAZ00001379873_PCSI_wgs_bam_PCSI_0102_Pa_P_ref.bam.c4gh
# EGAD00001004551/PART_2/EGAZ00001379874_PCSI_wgs_bam_PCSI_0102_Pa_P_tail.bam.c4gh
# EGAD00001004551/PART_1/EGAZ00001312063_COMPASS_wgs_PCSI_0664_Lv_M_526.bam.c4gh
# EGAD00001004551/PART_1/EGAZ00001312061_COMPASS_wgs_PCSI_0663_Lv_M_526.bam.c4gh
# EGAD00001004551/PART_1/EGAZ00001312042_COMPASS_wgs_PCSI_0646_Lv_M_526.bam.c4gh
# EGAD00001004551/PART_1/EGAZ00001312040_COMPASS_wgs_PCSI_0644_Lv_M_526.bam.c4gh
# EGAD00001004551/PART_1/EGAZ00001312028_COMPASS_wgs_PCSI_0632_Lv_M_526.bam.c4gh
# EGAD00001004551/PART_2/EGAZ00001379872_PCSI_wgs_bam_PCSI_0102_Pa_P_head.bam.c4gh
# EGAD00001003585/PART_1/EGAZ00001312063_COMPASS_wgs_PCSI_0664_Lv_M_526.bam.c4gh
# EGAD00001003585/PART_1/EGAZ00001312061_COMPASS_wgs_PCSI_0663_Lv_M_526.bam.c4gh
# EGAD00001003585/PART_1/EGAZ00001312042_COMPASS_wgs_PCSI_0646_Lv_M_526.bam.c4gh
# EGAD00001003585/PART_1/EGAZ00001312040_COMPASS_wgs_PCSI_0644_Lv_M_526.bam.c4gh
# EGAD00001003585/PART_1/EGAZ00001312028_COMPASS_wgs_PCSI_0632_Lv_M_526.bam.c4gh

# read report from EGA server sent by Niklas
dt.egaReport <- fread('../../EGAreport.csv')
dt.egaReport[, SampleID := gsub('(_unmapped_R[12])?.(bam|fastq)(.gz)?', '', FileName)]
dt.egaReport[, FileEnding := gsub('[A-Z]{4}[A-Z]?_[0-9]{4}_[A-Z][a-z]_[A-Z]_?([0-9]+)?([a-z]+)?(_unmapped_R[12])?', '', FileName)]
# status is 1 for all samples
dt.egaReport[, Status := NULL]

# merge dt.info table with ega report table
dt.merge <-merge(dt.info, dt.egaReport[FileEnding == '.bam'],
                 by.x = c("FileID", "FileName", "SampleID", "Bytes", "MD5"),
                 by.y = c("EGAF", "FileName", "SampleID", "Bytes", "CheckSum"),
                 all = T) %>% unique()
# yields table with 1219 rows (258 samples not from dt.info)

# read ICGC table sent by Sebastian
dt.icgc <- read_xlsx('../../ICGC_PACA_CA_WGS_849Bam_4Studies_SM.xlsx', 'All') %>% as.data.table()

# merge ICGC table with previous merge
# make Bytes numeric
dt.merge[, Bytes := as.numeric(Bytes)]
dt.merge2 <- merge(dt.merge, dt.icgc[, .(Study, EGA_Dataset, File_ID, File_name, CheckSum)],  #, Bytes, donor_id
      by.x = c("FileID", "FileName", "MD5"), #, "Bytes" 
      by.y = c("File_ID", "File_name", "CheckSum"), # , "Bytes"
      all = T) %>% unique()

# add batch = none instead of NA or no
dt.merge2[is.na(Batch) | Batch == 'no', Batch := 'none']

# some IDs are double because of the Bytes values -> do not use Bytes to merge?
# double_EGAF <- dt.merge2[, .N, by = FileID][N > 1, FileID]
# dt.merge2[FileID %in% double_EGAF]

# merge with dt.aspera 
dt.master <- merge(dt.merge2, dt.aspera, by = c("FileName"), all = T) %>% unique()
# contains 111 double entries with two uploads to Aspera
double_EGAF_master <- dt.master[, .N, by = FileID][N > 1, FileID]
dt.master[FileID %in% double_EGAF_master]

# TODO: replace NA values in Date, AnalysisFolverX, StatusLRZ, IgnoreFile, Study, EGAD
# TODO: resolve 103 entries with unmatching EGAD ids
dt.master[EGA_Dataset != EGAD]