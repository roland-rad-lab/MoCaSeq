#!/bin/bash
#SBATCH -J bam_chr_rename
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mail-type=end
#SBATCH --mem=12gb
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=marcus.wagner@tum.de
#SBATCH --export=NONE
#SBATCH --time=24:00:00
export OMP_NUM_THREADS=8

# script from Niklas to change chr names in bams using samtools

# get bam file
rawBAM=$1

# assert bam file suffix
if [[ $rawBAM != *.bam ]]
then
  echo "Wrong file ending. Expecting .bam, given $1"
  exit 1
fi

tmpBAM="${rawBAM/.bam/.tmp.bam}"
BAM="${rawBAM/.bam/.chrRenamed.bam}"

samtools view --threads $OMP_NUM_THREADS -o ${tmpBAM} ${rawBAM} `seq -f 'chr%g' 1 22` chrX chrY
samtools view -H ${tmpBAM} > header
sed '/random/d;/chrUn/d;/chrM/d;/chrEBV/d;/GL/d;/JH/d' -i header # remove all other chroms
sed -e 's/SN:chr/SN:/' -i header
samtools reheader header -i ${tmpBAM} > ${BAM}
samtools index ${BAM}
rm header ${tmpBAM}
