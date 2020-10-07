# INSTALL
working_directory=/mnt/3TBVol/MoCaSeq_WD/
mkdir -p ${working_directory} \
&& cd ${working_directory}

sudo curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh > script.deb.sh
sudo bash script.deb.sh
sudo apt-get install git-lfs
sudo git lfs install

sudo docker pull rolandradlab/mocaseq




# DEFINE
function set-title() {
  if [[ -z "$ORIG" ]]; then
    ORIG=$PS1
  fi
  TITLE="\[\e]2;$*\a\]"
  PS1=${ORIG}${TITLE}
}



# with custom script
# /home/rad/packages/MoCaSeq should be develop branch of github
sudo docker run \
-it --entrypoint=/bin/bash \
-v /mnt/3TBVol/MoCaSeq_WD:/var/pipeline/ \
-v /mnt/3TBVol/MoCaSeq_WD/MoCaSeq_ref:/var/pipeline/ref \
-v /home/rad/packages/MoCaSeq:/opt/MoCaSeq \
rolandradlab/mocaseq


sudo docker run \
-it --entrypoint=/bin/bash \
-v /mnt/3TBVol/MoCaSeq_WD:/var/pipeline/ \
-v /mnt/3TBVol/MoCaSeq_WD/ref_human_runtimeTest:/var/pipeline/ref \
-v /home/rad/packages/MoCaSeq:/opt/MoCaSeq \
rolandradlab/mocaseq


# for mouse with default script
sudo docker run \
-it --entrypoint=/bin/bash \
-v /mnt/3TBVol/MoCaSeq_WD:/var/pipeline/ \
-v /mnt/3TBVol/MoCaSeq_WD/MoCaSeq_ref_humanDNA_fromWS3:/var/pipeline/ref \
rolandradlab/mocaseq


sudo docker run \
-it --entrypoint=/bin/bash \
-v /mnt/3TBVol/MoCaSeq_WD:/var/pipeline/ \
-v /mnt/3TBVol/MoCaSeq_WD/MoCaSeq_ref:/var/pipeline/ref \
-v /home/rad/packages/MoCaSeq:/opt/MoCaSeq \
-v /mnt/3TBVol/data/AGRad_hPDAC_ICGC-PACA-CA-WGS:/var/fastqs/ \
rolandradlab/mocaseq






# DOCKER BUILD NIKLAS
sudo docker run \
-it --entrypoint=/bin/bash \
-v /media/rad/RAD8-2000/MoCaSeq_ICGC:/var/pipeline/ \
-v /mnt/3TBVol/MoCaSeq_WD/MoCaSeq_ref:/var/pipeline/ref \
-v /home/rad/packages/MoCaSeq:/opt/MoCaSeq \
-v /mnt/3TBVol/data/AGRad_hPDAC_ICGC-PACA-CA-WGS:/var/fastqs/ \
mocaseq-human

#cp: cannot stat 'temp/AGRad_hPDAC_ICGC-PACA-CA-WGS.Normal.cleaned.bam': No such file or directory
#cp: error reading 'temp/sambamba-pid4716-wuiz/AGRad_hPDAC_ICGC-PACA-CA-WGS.Normal.cleaned.bam.11': Input/output error



name=AGRad_hPDAC_ICGC-PACA-CA-WGS
# fastq_normal_1=/var/fastqs/EGAF00001709814/PCSI_0612_Si_R1.fastq.gz
# fastq_normal_2=/var/fastqs/EGAF00001709814/PCSI_0612_Si_R2.fastq.gz
# fastq_tumor_1=/var/fastqs/EGAF00001721862/PCSI_0612_Ag_M_526_R1.fastq.gz
# fastq_tumor_2=/var/fastqs/EGAF00001721862/PCSI_0612_Ag_M_526_R2.fastq.gz
fastq_normal_1=/var/fastqs/AGRad_hPDAC_ICGC-PACA-CA-WGS.Normal.R1.fastq.gz
fastq_normal_2=/var/fastqs/AGRad_hPDAC_ICGC-PACA-CA-WGS.Normal.R2.fastq.gz
fastq_tumor_1=/var/fastqs/AGRad_hPDAC_ICGC-PACA-CA-WGS.Tumor.R1.fastq.gz
fastq_tumor_2=/var/fastqs/AGRad_hPDAC_ICGC-PACA-CA-WGS.Tumor.R2.fastq.gz
bam_normal=
bam_tumor=
repeat_mapping=no

# name=AGRad_hPDAC_ICGC-PACA-CA-WGS-fromBAM
# fastq_normal_1=
# fastq_normal_2=
# fastq_tumor_1=
# fastq_tumor_2=
# bam_normal=/var/fastqs/EGAF00001709814/PCSI_0612_Si.bam
# bam_tumor=/var/fastqs/EGAF00001721862/PCSI_0612_Ag_M_526.bam
#repeat_mapping=yes

sequencing_type=WGS
quality_control=yes
threads=8
RAM=32
temp_dir=/var/pipeline/temp
artefact_type=none
filtering=hard
phred=
Mutect2=yes
Delly=no
runmode=MS
GATK=4.1.7.0
test=no
memstats=0
config_file=
species=Human
Titan=yes
Absolute=yes
Facets=yes
BubbleTree=yes





sudo docker run \
-it --entrypoint=/bin/bash \
-v /mnt/3TBVol/MoCaSeq_WD:/var/pipeline/ \
-v /mnt/3TBVol/MoCaSeq_WD/MoCaSeq_ref:/var/pipeline/ref \
-v /home/rad/packages/MoCaSeq:/opt/MoCaSeq \
-v /media/rad/RAD2-3000/AGRad/AGRad_hPDAC_ICGC-PACA-CA-WGS:/var/fastqs/ \
mocaseq_humanpipeline2

bash /opt/MoCaSeq/MoCaSeq.sh  \
  --name AGRad_hPDAC_ICGC-PACA-CA-WGS \
  --bam_tumor /var/fastqs/EGAF00001709814/PCSI_0612_Si_R.bam \
  --bam_normal /var/fastqs/EGAF00001721862/PCSI_0612_Ag_M_526.bam \
  --repeat_mapping yes \
  --sequencing_type WGS \
  --quality_control yes \
  --Titan yes \
  --Absolute yes \
  --Facets yes \
  --BubbleTree yes








# NO REMAPPING, USE BAMS AS THEY ARE
sudo docker run \
-it --entrypoint=/bin/bash \
-v /mnt/3TBVol/MoCaSeq_WD:/var/pipeline/ \
-v /mnt/3TBVol/MoCaSeq_WD/MoCaSeq_ref:/var/pipeline/ref \
-v /home/rad/packages/MoCaSeq:/opt/MoCaSeq \
-v /media/rad/RAD2-3000/AGRad/AGRad_hPDAC_ICGC-PACA-CA-WGS:/var/fastqs/ \
mocaseq_humanpipeline2

bash /opt/MoCaSeq/MoCaSeq.sh  \
  --name AGRad_hPDAC_ICGC-PACA-CA-WGS_noRemap \
  --bam_tumor /var/fastqs/EGAF00001709814/PCSI_0612_Si_R.bam \
  --bam_normal /var/fastqs/EGAF00001721862/PCSI_0612_Ag_M_526.bam \
  --repeat_mapping no \
  --sequencing_type WGS \
  --quality_control yes \
  --Titan yes \
  --Absolute yes \
  --Facets yes \
  --BubbleTree yes







sudo docker run \
-it --entrypoint=/bin/bash \
-v /mnt/3TBVol/MoCaSeq_WD:/var/pipeline/ \
-v /mnt/3TBVol/MoCaSeq_WD/MoCaSeq_ref:/var/pipeline/ref \
-v /home/rad/packages/MoCaSeq:/opt/MoCaSeq \
-v /mnt/3TBVol/data/AGKrackhardt_WES_forMoCaSeqTest:/var/fastqs/ \
mocaseq-human

bash /opt/MoCaSeq/MoCaSeq.sh  \
  --name AGKrackhardt_ZFQ9G4 \
  --fastq_normal_fw /var/fastqs/blood/AS-336552-LR-41660_R1.fastq.gz \
  --fastq_normal_rev /var/fastqs/blood/AS-336552-LR-41660_R2.fastq.gz \
  --fastq_tumor_fw /var/fastqs/tumor/AS-336553-LR-41660_R1.fastq.gz \
  --fastq_tumor_rev /var/fastqs/tumor/AS-336553-LR-41660_R2.fastq.gz \
  --species Human \
  --repeat_mapping yes \
  --sequencing_type WES \
  --quality_control yes \
  --Mutect2 yes \
  --Titan yes \
  --Absolute yes \
  --Facets yes \
  --BubbleTree yes







sudo docker run \
-it --entrypoint=/bin/bash \
-v /mnt/3TBVol/MoCaSeq_WD:/var/pipeline/ \
-v /mnt/3TBVol/MoCaSeq_WD/MoCaSeq_ref:/var/pipeline/ref \
-v /home/rad/packages/MoCaSeq:/opt/MoCaSeq \
-v /mnt/3TBVol/data/AGKrackhardt_WES_forMoCaSeqTest:/var/fastqs/ \
mocaseq_humanpipeline2


name=AGKrackhardt_ZFQ9G4
fastq_normal_1=/var/fastqs/blood/AS-336552-LR-41660_R1.fastq.gz
fastq_normal_2=/var/fastqs/blood/AS-336552-LR-41660_R2.fastq.gz
fastq_tumor_1=/var/fastqs/tumor/AS-336553-LR-41660_R1.fastq.gz
fastq_tumor_2=/var/fastqs/tumor/AS-336553-LR-41660_R2.fastq.gz
bam_normal=
bam_tumor=
repeat_mapping=yes
sequencing_type=WES
quality_control=no
threads=8
RAM=32
temp_dir=/var/pipeline/temp
artefact_type=none
filtering=hard
phred=
Mutect2=yes
Delly=no
runmode=MS
GATK=4.1.7.0
test=no
memstats=0
config_file=
species=Human
Titan=yes
Absoute=yes
Facets=yes
BubbleTree=yes
