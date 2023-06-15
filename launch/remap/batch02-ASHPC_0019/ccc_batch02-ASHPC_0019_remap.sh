#!/bin/bash
#SBATCH -J ASHPC_0019-remap
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --clusters=mpp3
#SBATCH --partition=mpp3_batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=30
#SBATCH --get-user-env
#SBATCH --mail-type=end
#SBATCH --mail-user=marcus.wagner@tum.de
#SBATCH --export=NONE
#SBATCH --time=48:00:00

# SLURM script to remap ASHPC_0019 samples for MoCaSeq from hg19 to hg38 using charliecloud container (aka ccc) on a single node of CoolMUC-3 cluster.
# @author: marcus.wagner@tum.de.
# Cluster options:
# - cm2_tiny: max 72h, 56GB, 28 CPUs, jobs(10/50), nextflow compatible 
# - serial_long: max 480h, 56GB, 28 CPUs, jobs(dynamic,250), nextflow compatible 
# - mpp3_batch: max 48h, 90GB, 64 CPUs, jobs(50,dynamic), only bash version of MoCaSeq

module load slurm_setup

# set container path
cccDir=${HOME}/images-live/mocaseq2

# specify mount paths
workingDir=/gpfs/scratch/pn29ya/$USER/mocaseq-slurm/remap
mocaseqDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/software/MoCaSeq/
bamDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS//input/GRCh37_bam/batch02/
referencesDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/references/bashMoCaSeq/
mkdir -p $workingDir
ln -s -t $workingDir $referencesDir

# specify sample
sample=ASHPC_0019_Pa_P_5262
bamDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS//input/GRCh37_bam/batch02/
bamName=EGAZ00001379792_PCSI_wgs_bam_ASHPC_0019_Pa_P_5262.bam
bamType="Tumor"
# submit subjob for sample remapping
srun --ntasks=1 --exclusive --mem 45568mb -J $sample -o ./%x.%j.%N.out ${mocaseqDir}/launch/ccc_remap_wrapper.sh -ccc $cccDir -wd $workingDir -m $mocaseqDir -bd $bamDir -bf $bamName -s $sample -rd $referencesDir -t $bamType > ${sample}-remap.out & 
sleep 4

# specify sample
sample=ASHPC_0019_Pa_R
bamDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS//input/GRCh37_bam/batch02/
bamName=EGAZ00001379793_PCSI_wgs_bam_ASHPC_0019_Pa_R.bam
# submit subjob for sample remapping
srun --ntasks=1 --exclusive --mem 45568mb -J $sample -o ./%x.%j.%N.out ch-run $cccDir --no-home --set-env=sample=${sample} -w --no-passwd \
--bind ${workingDir}:/var/pipeline/ \
--bind ${mocaseqDir}:/opt/MoCaSeq/ \
--bind ${bamDir}:/var/raw-bams/ \
--bind $referencesDir:$referencesDir \
-- \
/opt/MoCaSeq/MoCaSeq_LRZ_remap_postprocess.sh \
-nb /var/raw-bams/${bamName} \
--name ${sample} \
--species Human \
--repeat_mapping yes \
--sequencing_type WGS \
--quality_control yes \
--threads 30 \
--RAM 44 \
--GATKVersion 4.1.7.0 \
--filtering soft \
--artefact yes > ${sample}-remap.out & 

wait # for completion of background tasks

