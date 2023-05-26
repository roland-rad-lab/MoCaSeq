#!/bin/bash
# Test script to run MoCaSeq in charliecloud using single SLURM job

#SBATCH -J slurm-ccc-test
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --get-user-env
#SBATCH --mail-type=end
#SBATCH --mem=800mb
#SBATCH --mail-user=marcus.wagner@tum.de
#SBATCH --export=NONE
#SBATCH --time=08:00:00

module load slurm_setup
# module load ...
ccc_path=${HOME}/images-live/mocaseq2/

# here ch-run ... java -Xmx8G -jar /opt/gatk-4.2.0.0/gatk.jar Mutect2 --help
ch-run $ccc_path -- java -Xmx8G -jar /opt/gatk-4.2.0.0/gatk.jar Mutect2 --help > $HOME/mutect_help.txt


