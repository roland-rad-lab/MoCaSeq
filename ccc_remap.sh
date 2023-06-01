#!/bin/bash
# SLURM template to remap samples for MoCaSeq from hg19 to hg38 using charliecloud container (aka ccc)
# @author: marcus.wagner@tum.de

#SBATCH -J remapping-test
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

# set container path
# for shared container use /gpfs/scratch/pn29ya/ga89tog2/charliecloud-containers/mocaseq2
ccc_path=${HOME}/images-live/mocaseq2

# specify mount paths
working_directory=/gpfs/scratch/pn29ya/ga89tog2/mocaseq-slurm
# ref_directory=/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/reference_bash # copied reference dir to ${working_directory}/ref
script_directory=/dss/dsshome1/0F/ga89tog2/.nextflow/assets/roland-rad-lab/MoCaSeq
bam_path_prefix=/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/projects/hPDAC/ICGC_PACA_CA_WGS/ega_download/debug_bams/
# TODO out dir

# specify samples
sample=PCSI_0357_St_R_1percent

# MoCaSeq-remap call inside charliecloud container
ch-run $ccc_path --no-home --set-env=sample=${sample} -w --no-passwd \
--bind ${working_directory}:/var/pipeline/ \
--bind ${script_directory}:/opt/MoCaSeq/ \
--bind ${bam_path_prefix}:/var/raw-bams/ \
-- \
/opt/MoCaSeq/MoCaSeq_LRZ_remap.sh \
-nb /var/raw-bams/${sample}.bam \
--name ${sample} \
--species Human \
--repeat_mapping yes \
--sequencing_type WGS \
--quality_control no \
--threads 40 \
--RAM 40 \
--GATKVersion 4.1.7.0 \
--filtering soft \
-qc yes \
--artefact yes

