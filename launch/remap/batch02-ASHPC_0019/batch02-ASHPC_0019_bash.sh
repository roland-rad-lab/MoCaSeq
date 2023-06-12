#!/bin/bash
#SBATCH -J batch02-ASHPC_0019-remap
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --clusters=mpp3
#SBATCH --partition=mpp3_batch
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=32
#SBATCH --mem=89gb
#SBATCH --get-user-env
#SBATCH --mail-type=end
#SBATCH --mail-user=marcus.wagner@tum.de
#SBATCH --export=NONE

# SLURM script for remapping sample group batch02-ASHPC_0019 .bam files hg19 -> hg38.
# Execution on a single node of CoolMUC-3 cluster.
# ClusterInfo: mpp3_batch: max 48h, 90GB, 64 CPUs, jobs(50,dynamic)
# See https://doku.lrz.de/job-processing-on-the-linux-cluster-10745970.html
# !!! Please check:
# - you have the omcaseq2 charldiecloud container in $HOME/images-live
# - all paths are accessible, expecially the reference path 

module load slurm_setup

# set container path
cccDir=${HOME}/images-live/mocaseq2

# specify mount paths
workingDir=/gpfs/scratch/pn29ya/$USER/mocaseq-slurm/remap
# ref_directory=/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/reference_bash # copied reference dir to ${working_directory}/ref
mocaseqDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/software/MoCaSeq/
referencesDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/references/bashMoCaSeq/
# make sure the working dir exists and create a symbolic link to the references named "ref" inside working directory
mkdir -p $workingDir
mkdir -p /gpfs/scratch/pn29ya/$USER/tmp
ln -s $referencesDir ref
mv ref $workingDir

# specify sample
sample=ASHPC_0019_Pa_P_5262
bamDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS//input/GRCh37_bam/batch02
bamName=EGAZ00001379792_PCSI_wgs_bam_ASHPC_0019_Pa_P_5262.bam
bamType="Tumor"
# submit subjob for sample remapping
srun --ntasks=1 --exclusive --mem 44.5gb -J $sample -o ./%x.%j.%N.out ${mocaseqDir}/launch/ccc_remap_wrapper.sh -ccc $cccDir -wd $workingDir -m $mocaseqDir -bd $bamDir -bf $bamName -s $sample -rd $referencesDir -t $bamType > ${sample}-remap.out & 
sleep 4

# specify sample
sample=ASHPC_0019_Pa_R
bamDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS//input/GRCh37_bam/batch02
bamName=EGAZ00001379793_PCSI_wgs_bam_ASHPC_0019_Pa_R.bam
bamType="Normal"
# submit subjob for sample remapping
srun --ntasks=1 --exclusive --mem 44.5gb -J $sample -o ./%x.%j.%N.out ${mocaseqDir}/launch/ccc_remap_wrapper.sh -ccc $cccDir -wd $workingDir -m $mocaseqDir -bd $bamDir -bf $bamName -s $sample -rd $referencesDir -t $bamType > ${sample}-remap.out & 
sleep 4

wait # for completion of background tasks
