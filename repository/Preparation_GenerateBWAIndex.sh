#!/bin/bash

##########################################################################################
##
## Preparation_GenerateBWAIndex.sh
##
## Uses a downloaded .fasta-file to generate the BWA index needed for contig-aware mapping.
##
##########################################################################################

VersionMouse=$1
fasta_to_fastq=$2

sed -ri "s/>CM[0-9\.]* Mus musculus chromosome ([0-9XY]*).*/>\1/g" Genomes/$VersionMouse/$VersionMouse.fna
samtools faidx Genomes/$VersionMouse/$VersionMouse.fna
samtools faidx Genomes/$VersionMouse/$VersionMouse.fna $(grep -P -o '^>[A-Z].*\.\d*' Genomes/$VersionMouse/$VersionMouse.fna | sed 's/>//g') > Genomes/$VersionMouse/"$VersionMouse"_alt.fna
samtools faidx Genomes/$VersionMouse/$VersionMouse.fna $(grep -P -o '>[0-9XY]+' Genomes/$VersionMouse/$VersionMouse.fna | sed 's/>//g') > Genomes/$VersionMouse/"$VersionMouse"_primary.fna

perl $fasta_to_fastq/fasta_to_fastq.pl Genomes/$VersionMouse/"$VersionMouse"_alt.fna > Genomes/$VersionMouse/haplotypes.fastq

mkdir Genomes/$VersionMouse/primary_index

bwa index -p Genomes/$VersionMouse/primary_index/bwa_index -a bwtsw Genomes/$VersionMouse/"$VersionMouse"_primary.fna

bwa mem -t $threads Genomes/$VersionMouse/primary_index/bwa_index Genomes/$VersionMouse/haplotypes.fastq > Genomes/$VersionMouse/alt_mapping_unsorted.sam
samtools sort -@ $threads -O sam -o Genomes/$VersionMouse/alt_mapping_nosup.sam Genomes/$VersionMouse/alt_mapping_unsorted.sam
samtools view -h -o Genomes/$VersionMouse/alt_mapping.sam -F 0x800 Genomes/$VersionMouse/alt_mapping_nosup.sam

mkdir Genomes/$VersionMouse/tmp

grep '^@SQ' Genomes/$VersionMouse/alt_mapping.sam > Genomes/$VersionMouse/tmp/header
grep '^[^@]' Genomes/$VersionMouse/alt_mapping.sam > Genomes/$VersionMouse/tmp/data
cut -f 1 Genomes/$VersionMouse/tmp/data > Genomes/$VersionMouse/tmp/name
cut -f 3-4 Genomes/$VersionMouse/tmp/data > Genomes/$VersionMouse/tmp/ref_pos
cut -f 6 Genomes/$VersionMouse/tmp/data > Genomes/$VersionMouse/tmp/cigar
grep -P -o 'NM:i:\d*' Genomes/$VersionMouse/tmp/data > Genomes/$VersionMouse/tmp/dist

length=$(grep -c '^[^@]' Genomes/$VersionMouse/alt_mapping.sam)

for i in $(seq 1 $length);
do
	echo '0' >> Genomes/$VersionMouse/tmp/zero
	echo '*' >> Genomes/$VersionMouse/tmp/asterisk
	echo '255' >> Genomes/$VersionMouse/tmp/mapq
done

mkdir Genomes/$VersionMouse/bwa_index

paste Genomes/$VersionMouse/tmp/name Genomes/$VersionMouse/tmp/zero Genomes/$VersionMouse/tmp/ref_pos Genomes/$VersionMouse/tmp/mapq Genomes/$VersionMouse/tmp/cigar Genomes/$VersionMouse/tmp/asterisk Genomes/$VersionMouse/tmp/zero Genomes/$VersionMouse/tmp/zero Genomes/$VersionMouse/tmp/asterisk Genomes/$VersionMouse/tmp/asterisk Genomes/$VersionMouse/tmp/dist | cat Genomes/$VersionMouse/tmp/header - > Genomes/$VersionMouse/bwa_index/$VersionMouse'.alt'

bwa index -p Genomes/$VersionMouse/bwa_index/$VersionMouse -a bwtsw Genomes/$VersionMouse/$VersionMouse.fna

rm Genomes/$VersionMouse/alt_mapping.sam Genomes/$VersionMouse/alt_mapping_nosup.sam Genomes/$VersionMouse/alt_mapping_unsorted.sam
rm Genomes/$VersionMouse/"$VersionMouse"_alt.fna
rm Genomes/$VersionMouse/"$VersionMouse"_primary.fna
rm Genomes/$VersionMouse/haplotypes.fastq
rm -r Genomes/$VersionMouse/primary_index
rm -r Genomes/$VersionMouse/tmp