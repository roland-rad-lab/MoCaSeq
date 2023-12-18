#!/usr/bin/sh

# VCF to MAF file conversion in charliecloud container
# generate MAF file with variant effect predicted (VEP) (i.e. frameshift -> Frame_Shift_Del)
# sh $repository_dir/SNV_RunVEP.sh $name $config_file $species Mutect2 $runmode $types

ccc_path=${HOME}/images-live/mocaseq2
working_directory=/gpfs/scratch/pn29ya/ga89tog2/mocaseq-slurm
script_directory=/dss/dsshome1/0F/ga89tog2/.nextflow/assets/roland-rad-lab/MoCaSeq
output_directory=/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/projects/hPDAC/ICGC_PACA_CA_WGS/output/batch04

# TODO detect file names of all matched.Mutect2.txt vcf files and apply SNV_RunVEP.sh for all of them.
files=$(find output_directory -name *Mutect2.txt)
config_file=/opt/MoCaSeq/config.sh
species=Human
runmode=matched
types=Tumor

for name in files;
do
# MoCaSeq-remap call inside charliecloud container
# TODO PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/bin
ch-run $ccc_path --no-home -w --no-passwd \
--bind ${working_directory}:/var/pipeline/ \
--bind ${script_directory}:/opt/MoCaSeq/ \
--bind ${output_directory}:/var/output/ \
--cd /var/output/ \
-- \
/opt/MoCaSeq/repository/SNV_RunVEP.sh $name $config_file $species Mutect2 $runmode $types

done

