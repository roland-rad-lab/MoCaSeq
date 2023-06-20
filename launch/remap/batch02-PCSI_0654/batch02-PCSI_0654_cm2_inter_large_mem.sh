#!/bin/bash
#SBATCH -J batch02-PCSI_0654-remap
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --clusters=inter
#SBATCH --partition=cm2_inter_large_mem
#SBATCH --nodes=1
#SBATCH --time=96:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=121856mb
#SBATCH --get-user-env
#SBATCH --mail-type=end
#SBATCH --mail-user=marcus.wagner@tum.de
#SBATCH --export=NONE

# SLURM script for remapping sample group batch02-PCSI_0654 .bam files hg19 -> hg38.
# Execution on a single node of CoolMUC-2 cluster cm2_inter_large_mem partition.
# ClusterInfo: cm2_inter_large_mem: max 96h, 119GB, 96 CPUs, jobs(1,2)
# See https://doku.lrz.de/job-processing-on-the-linux-cluster-10745970.html
# !!! Please check:
# - you have the mocaseq2 charldiecloud container in $HOME/images-live
# - all paths are accessible, expecially the reference path 

module load slurm_setup

# set container path
cccDir=${HOME}/images-live/mocaseq2

# specify mount paths
workingDir=/gpfs/scratch/pn29ya/$USER/mocaseq-slurm/remap
mocaseqDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/software/MoCaSeq/
referencesDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/references/bashMoCaSeq/
# make sure the working dir exists and create a symbolic link to the references named "ref" inside working directory
mkdir -p $workingDir
mkdir -p /gpfs/scratch/pn29ya/$USER/tmp
ln -s $referencesDir ref
mv ref $workingDir

# specify sample
sample=PCSI_0654_Lv_M_526
bamDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS//input/GRCh37_bam/batch02
bamName=EGAZ00001312052_COMPASS_wgs_PCSI_0654_Lv_M_526.bam
bamType="Tumor"
# submit subjob for sample remapping
srun --ntasks=1 --exclusive --mem 60416mb -J $sample -o ./%x.%j.%N.out ${mocaseqDir}/launch/ccc_remap_wrapper.sh -ccc $cccDir -wd $workingDir -m $mocaseqDir -bd $bamDir -bf $bamName -s $sample -rd $referencesDir -t $bamType -r 59 -@ 48 > ${sample}-remap.out & 
sleep 4

# specify sample
sample=PCSI_0654_Lv_M_5262
bamDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS//input/GRCh37_bam/batch02
bamName=EGAZ00001383193_PCSI_wgs_bam_PCSI_0654_Lv_M_5262.bam
bamType="Tumor"
# submit subjob for sample remapping
srun --ntasks=1 --exclusive --mem 60416mb -J $sample -o ./%x.%j.%N.out ${mocaseqDir}/launch/ccc_remap_wrapper.sh -ccc $cccDir -wd $workingDir -m $mocaseqDir -bd $bamDir -bf $bamName -s $sample -rd $referencesDir -t $bamType -r 59 -@ 48 > ${sample}-remap.out & 
sleep 4

wait # for completion of background tasks