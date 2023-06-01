#!/bin/bash
# SLURM template for MoCaSeq in charliecloud container (aka ccc)
# @author: marcus.wagner@tum.de

#SBATCH -J analysis-test
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --clusters=cm2
#SBATCH --partition=cm2_std
#SBATCH --get-user-env
#SBATCH --mail-type=end
#SBATCH --mem=40gb
#SBATCH --mail-user=marcus.wagner@tum.de
#SBATCH --export=NONE
#SBATCH --time=08:00:00

module load slurm_setup

# set container path
# for shared container use /gpfs/scratch/pn29ya/ga89tog2/charliecloud-containers/mocaseq2
ccc_path=${HOME}/images-live/mocaseq2

# specify mount path
working_directory=/gpfs/scratch/pn29ya/ga89tog2/mocaseq-slurm
# ref_directory=/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/reference_bash # copied reference dir to ${working_directory}/ref
script_directory=/dss/dsshome1/0F/ga89tog2/.nextflow/assets/roland-rad-lab/MoCaSeq
# /dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/projects/hPDAC/ICGC_PACA_CA_WGS/input/GRCh38.p12_bam/batch01/
bam_path_prefix=/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/projects/hPDAC/ICGC_PACA_CA_WGS/ega_download/debug_bams/bam_remap
# TODO out dir

# specify samples
sampleT=PCSI_0357_Pa_P_5262_1percent
sampleN=PCSI_0357_St_R_1percent


# MoCaSeq call inside charliecloud container
# TODO: adapt MoCaSeq paths & parameters to sample call
ch-run $ccc_path --no-home --set-env=sampleT=${sampleT} --set-env=sampleT=${sampleN} -w --no-passwd \
--bind ${working_directory}:/var/pipeline/ \
--bind ${script_directory}:/opt/MoCaSeq/ \
--bind ${bam_path_prefix}:/var/raw-bams/ \
-- \
/opt/MoCaSeq/MoCaSeq_COMPASS_WGS.sh \
-tb /var/raw-bams/${sampleT}/${sampleT}.Tumor.bam \
-nb /var/raw-bams/${sampleN}/${sampleN}.Normal.bam \
--name ${sampleT} \
--species Human \
--repeat_mapping no \
--sequencing_type WGS \
--quality_control no \
--threads 40 \
--RAM 40 \
--GATKVersion 4.1.7.0 \
--filtering soft \
--artefact yes \
--Mutect2 yes \
--CNVKit yes \
--Delly no \
--BubbleTree yes \
--Absolute yes \
--Facets yes \
--Titan yes
# --para
