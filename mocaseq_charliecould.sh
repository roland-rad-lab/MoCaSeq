#!/bin/bash
# Test script to run MoCaSeq in charliecloud using single SLURM job

#SBATCH -J RAMP_0008-test
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --clusters=cm2_tiny
#SBATCH --get-user-env
#SBATCH --mail-type=end
#SBATCH --mem=40gb
#SBATCH --mail-user=marcus.wagner@tum.de
#SBATCH --export=NONE
#SBATCH --time=08:00:00

module load slurm_setup
# module load ...

# modify SLURM cluster to run on
export SLURM\_CLUSTERS="cm2"

# set container path
ccc_path=${HOME}/images-live/mocaseq2/

# specify mount path
working_directory=/gpfs/scratch/pn29ya/ga89tog2/mocaseq-slurm
ref_directory=/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/reference
script_directory=/dss/dsshome1/0F/ga89tog2/.nextflow/assets/roland-rad-lab/MoCaSeq
# TODO: asseble full bam path in pipeline: RAMP_0008_Pa_P/results/bam_remap/
bam_path_prefix=/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/projects/hPDAC/ICGC_PACA_CA_WGS/input/GRCh38.p12_bam/batch01/
# TODO out dir

# specify samples
sampleT=RAMP_0008_Pa_P
sampleN=RAMP_0008_Mu_R

# MoCaSeq call inside charliecloud container
ch-run $ccc_path --no-home --set-env -w --no-passwd \
--bind ${working_directory}:/var/wd/ \
--bind ${working_directory}/tmp:/var/pipeline/temp/ \
--bind ${ref_directory}:/var/pipeline/ref/ \
--bind ${script_directory}:/opt/MoCaSeq/ \
--bind ${bam_path_prefix}:/var/pipeline/raw/ \
-- \
/opt/MoCaSeq/MoCaSeq_LRZ.sh \
-tb /var/pipeline/raw/${sampleT}/results/bam_remap/${sampleT}.Tumor.bam \
-nb /var/pipeline/raw/${sampleN}/results/bam_remap/${sampleN}.Normal.bam \
--name ${sampleT} \
--species Human \
--repeat_mapping no \
--sequencing_type WGS \
--quality_control no \
--threads 28 \
--RAM 40 \
--filtering soft \
--artefact yes \
--Mutect2 yes \
--CNVKit yes \
--para