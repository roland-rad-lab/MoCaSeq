#!/bin/bash
#SBATCH -J PCSI_0705-remap
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --reservation=gen_seq
#SBATCH --mem=48gb
#SBATCH --get-user-env
#SBATCH --mail-type=end
#SBATCH --mail-user=marcus.wagner@tum.de
#SBATCH --export=NONE
#SBATCH --time=96:00:00

# SLURM template for running MoCaSeq nextflow pipeline on a single node of CoolMUC-3 cluster.
# @author: marcus.wagner@tum.de
# For this approach, do NOT use the slurm nextflow profile to have all processes executed "locally" on the allocated cluster node.
# Cluster options:
# - cm2_tiny: max 72h, 56GB, 28 CPUs, jobs(10/50), nextflow compatible 
# - serial_long: max 480h, 56GB, 28 CPUs, jobs(dynamic,250), nextflow compatible 
# - mpp3_batch: max 48h, 90GB, 64 CPUs, jobs(50,dynamic), only bash version of MoCaSeq
# 
# Important:
# - make sure to have $SLURM_CLUSTERS set to the cluster you want to use!
# - make sure java, nextflow and charliecloud are accessible in $PATH, $JAVA_HOME and $JAVA_CMD
# - make sure you have the charldiecloud containers in $HOME/images-live

module load slurm_setup

# specify output/project dir
projectDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS

# define MoCaSeq version
mocaseqDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/software/MoCaSeq

# specify working dir
workingDir=/gpfs/scratch/pn29ya/$USER/mocaseq-nextflow/remap
mkdir -p $workingDir

# nextflow call with charliecloud profile 
nextflow run ${mocaseqDir}/main.nf -profile charliecloud -entry MAP -work-dir ${workingDir} --output_base ${projectDir}/input --genome_build.human GRCh38.p12 --custom_config_base ${mocaseqDir}/conf --custom_config_version mocaseq-lrz --input ${mocaseqDir}/input/remap/batch02-PCSI_0705.tsv -resume

