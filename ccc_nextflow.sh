#!/bin/bash
# SLURM template for MoCaSeq in charliecloud container (aka ccc)
# @author: marcus.wagner@tum.de

#SBATCH -J analysis-test
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --clusters=inter
#SBATCH --partition=cm2_inter
#SBATCH --nodes=1
#SBATCH --get-user-env
#SBATCH --mail-type=end
#SBATCH --mem=40gb
#SBATCH --mail-user=marcus.wagner@tum.de
#SBATCH --export=NONE
#SBATCH --time=02:00:00

module load slurm_setup

# specify mount path
working_directory=/gpfs/scratch/pn29ya/ga89tog2/ga89tog2/test/

# nextflow call with charliecloud profile 
nextflow run main_new.nf -profile charliecloud -work-dir ${working_directory} --custom_config_base ./conf --custom_config_version mocaseq-lrz


