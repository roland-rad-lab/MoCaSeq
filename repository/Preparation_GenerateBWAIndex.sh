#!/bin/bash

##########################################################################################
##
## Preparation_GenerateBWAIndex.sh
##
## Uses a downloaded .fasta-file to generate the BWA index needed for contig-aware mapping.
##
##########################################################################################

VersionMouse=$1
config_file=$2

species=Mouse

. $config_file

sed -ri "s/>CM[0-9\.]* Mus musculus chromosome ([0-9XY]*).*/>\1/g" ref/$VersionMouse/$VersionMouse.fna
samtools faidx ref/$VersionMouse/$VersionMouse.fna
samtools faidx ref/$VersionMouse/$VersionMouse.fna $(grep -P -o '^>[A-Z].*\.\d*' ref/$VersionMouse/$VersionMouse.fna | sed 's/>//g') > ref/$VersionMouse"/"$VersionMouse"_alt.fna"
samtools faidx ref/$VersionMouse/$VersionMouse.fna $(grep -P -o '>[0-9XY]+' ref/$VersionMouse/$VersionMouse.fna | sed 's/>//g') > ref/$VersionMouse"/"$VersionMouse"_primary.fna"

perl $fasta_to_fastq"/fasta_to_fastq.pl" "ref/"$VersionMouse"/"$VersionMouse"_alt.fna" > "ref/"$VersionMouse"/haplotypes.fastq"

mkdir ref/$VersionMouse/primary_index

bwa index -p "ref/"$VersionMouse"/primary_index/bwa_index" -a bwtsw "ref/"$VersionMouse"/"$VersionMouse"_primary.fna"

bwa mem ref/$VersionMouse/primary_index/bwa_index ref/$VersionMouse/haplotypes.fastq > ref/$VersionMouse/alt_mapping_unsorted.sam

samtools sort -O sam -o ref/$VersionMouse/alt_mapping_nosup.sam ref/$VersionMouse/alt_mapping_unsorted.sam

samtools view -h -o ref/$VersionMouse/alt_mapping.sam -F 0x800 ref/$VersionMouse/alt_mapping_nosup.sam

mkdir ref/$VersionMouse/tmp

grep '^@SQ' ref/$VersionMouse/alt_mapping.sam > ref/$VersionMouse/tmp/header
grep '^[^@]' ref/$VersionMouse/alt_mapping.sam > ref/$VersionMouse/tmp/data
cut -f 1 ref/$VersionMouse/tmp/data > ref/$VersionMouse/tmp/name
cut -f 3-4 ref/$VersionMouse/tmp/data > ref/$VersionMouse/tmp/ref_pos
cut -f 6 ref/$VersionMouse/tmp/data > ref/$VersionMouse/tmp/cigar
grep -P -o 'NM:i:\d*' ref/$VersionMouse/tmp/data > ref/$VersionMouse/tmp/dist

length=$(grep -c '^[^@]' ref/$VersionMouse/alt_mapping.sam)

for i in $(seq 1 $length);
do
	echo '0' >> ref/$VersionMouse/tmp/zero
	echo '*' >> ref/$VersionMouse/tmp/asterisk
	echo '255' >> ref/$VersionMouse/tmp/mapq
done

mkdir ref/$VersionMouse/bwa_index

paste ref/$VersionMouse/tmp/name ref/$VersionMouse/tmp/zero ref/$VersionMouse/tmp/ref_pos ref/$VersionMouse/tmp/mapq ref/$VersionMouse/tmp/cigar ref/$VersionMouse/tmp/asterisk ref/$VersionMouse/tmp/zero ref/$VersionMouse/tmp/zero ref/$VersionMouse/tmp/asterisk ref/$VersionMouse/tmp/asterisk ref/$VersionMouse/tmp/dist | cat ref/$VersionMouse/tmp/header - > ref/$VersionMouse/bwa_index/$VersionMouse'.alt'

bwa index -p ref/$VersionMouse/bwa_index/$VersionMouse -a bwtsw ref/$VersionMouse/$VersionMouse.fna

rm ref/$VersionMouse/alt_mapping.sam ref/$VersionMouse/alt_mapping_nosup.sam ref/$VersionMouse/alt_mapping_unsorted.sam
rm "ref/"$VersionMouse"/"$VersionMouse"_alt.fna"
rm "ref/"$VersionMouse"/"$VersionMouse"_primary.fna"
rm ref/$VersionMouse/haplotypes.fastq
rm -r ref/$VersionMouse/primary_index
rm -r ref/$VersionMouse/tmp