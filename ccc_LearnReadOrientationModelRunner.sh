#!/bin/bash
#SBATCH -o %x.%j.%N.out
#SBATCH -e %x.%j.%N.err
#SBATCH -D ./
#SBATCH -J LearnReadModel
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --clusters=mpp3
#SBATCH --partition=mpp3_batch

module load slurm_setup

# run charliecloud container to do Mutect2 post processing
# set container path
cccDir=${HOME}/images-live/mocaseq2

# specify mount paths
workingDir=/gpfs/scratch/pn29ya/$USER/mocaseq-slurm/mocaseq/
mocaseqDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/software/MoCaSeq/
outDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/output/GRCh38.p12/
referencesDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/references/bashMoCaSeq/
# make sure the working dir exists and create a symbolic link to the references named "ref" inside working directory
mkdir -p $workingDir
mkdir -p /gpfs/scratch/pn29ya/$USER/tmp
ln -s $referencesDir ref
mv ref $workingDir

# postprocessing variables
name=$1
type=$2
startDir=$3
echo "${name} Mutect2 learn read model for $type case"

# MoCaSeq call inside charliecloud container
ch-run $cccDir --no-home --set-env=name=${name} -w --no-passwd \
--bind ${workingDir}:/var/pipeline/ \
--bind ${mocaseqDir}:/opt/MoCaSeq/ \
--bind ${outDir} \
--bind ${referencesDir} \
-c $startDir \
-- /bin/bash /opt/MoCaSeq/launch/ccc_LearnReadOrientationModel.sh -t $type
