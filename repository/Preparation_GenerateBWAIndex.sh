#!/bin/bash

##########################################################################################
##
## Preparation_GenerateBWAIndex.sh
##
## Uses a downloaded .fasta-file to generate the BWA index needed for contig-aware mapping.
##
##########################################################################################

VersionGenome=$1
config_file=$2
species=$3

. $config_file

# change header names for each species
if [ $species = 'Mouse' ]; then
	sed -ri "s/>CM[0-9\.]* Mus musculus chromosome ([0-9XY]*).*/>\1/g" ref/$VersionGenome/$VersionGenome.fna
elif [ $species = 'Human' ]; then
	sed -ri "s/>CM[0-9\.]* Homo sapiens chromosome ([0-9XY]*).*/>\1/g" ref/$VersionGenome/$VersionGenome.fna
else echo "Invalid species input (${species}). Choose Mouse or Human"; exit 1
fi

samtools faidx ref/${VersionGenome}/${VersionGenome}.fna
samtools faidx ref/${VersionGenome}/${VersionGenome}.fna $(grep -P -o '^>[A-Z].*\.\d*' ref/${VersionGenome}/${VersionGenome}.fna | sed 's/>//g') > ref/${VersionGenome}"/"${VersionGenome}"_alt.fna"
samtools faidx ref/${VersionGenome}/${VersionGenome}.fna $(grep -P -o '>[0-9XY]+' ref/${VersionGenome}/${VersionGenome}.fna | sed 's/>//g') > ref/${VersionGenome}"/"${VersionGenome}"_primary.fna"

perl /opt/bin/fasta_to_fastq.pl "ref/"${VersionGenome}"/"${VersionGenome}"_alt.fna" > "ref/"${VersionGenome}"/haplotypes.fastq"

mkdir -p ref/${VersionGenome}/primary_index

bwa index -p "ref/"${VersionGenome}"/primary_index/bwa_index" -a bwtsw "ref/"${VersionGenome}"/"${VersionGenome}"_primary.fna"

bwa mem ref/${VersionGenome}/primary_index/bwa_index ref/${VersionGenome}/haplotypes.fastq > ref/${VersionGenome}/alt_mapping_unsorted.sam

samtools sort -O sam -o ref/${VersionGenome}/alt_mapping_nosup.sam ref/${VersionGenome}/alt_mapping_unsorted.sam

samtools view -h -o ref/${VersionGenome}/alt_mapping.sam -F 0x800 ref/${VersionGenome}/alt_mapping_nosup.sam

mkdir -p ref/${VersionGenome}/tmp

grep '^@SQ' ref/${VersionGenome}/alt_mapping.sam > ref/${VersionGenome}/tmp/header
grep '^[^@]' ref/${VersionGenome}/alt_mapping.sam > ref/${VersionGenome}/tmp/data
cut -f 1 ref/${VersionGenome}/tmp/data > ref/${VersionGenome}/tmp/name
cut -f 3-4 ref/${VersionGenome}/tmp/data > ref/${VersionGenome}/tmp/ref_pos
cut -f 6 ref/${VersionGenome}/tmp/data > ref/${VersionGenome}/tmp/cigar
grep -P -o 'NM:i:\d*' ref/${VersionGenome}/tmp/data > ref/${VersionGenome}/tmp/dist

length=$(grep -c '^[^@]' ref/${VersionGenome}/alt_mapping.sam)

for i in $(seq 1 $length);
do
	echo '0' >> ref/${VersionGenome}/tmp/zero
	echo '*' >> ref/${VersionGenome}/tmp/asterisk
	echo '255' >> ref/${VersionGenome}/tmp/mapq
done

mkdir -p ref/${VersionGenome}/bwa_index

paste ref/${VersionGenome}/tmp/name ref/${VersionGenome}/tmp/zero ref/${VersionGenome}/tmp/ref_pos ref/${VersionGenome}/tmp/mapq ref/${VersionGenome}/tmp/cigar ref/${VersionGenome}/tmp/asterisk ref/${VersionGenome}/tmp/zero ref/${VersionGenome}/tmp/zero ref/${VersionGenome}/tmp/asterisk ref/${VersionGenome}/tmp/asterisk ref/${VersionGenome}/tmp/dist | cat ref/${VersionGenome}/tmp/header - > ref/${VersionGenome}/bwa_index/${VersionGenome}'.alt'

bwa index -p ref/${VersionGenome}/bwa_index/${VersionGenome} -a bwtsw ref/${VersionGenome}/${VersionGenome}.fna

rm ref/${VersionGenome}/alt_mapping.sam ref/${VersionGenome}/alt_mapping_nosup.sam ref/${VersionGenome}/alt_mapping_unsorted.sam
rm "ref/"${VersionGenome}"/"${VersionGenome}"_alt.fna"
rm "ref/"${VersionGenome}"/"${VersionGenome}"_primary.fna"
rm ref/${VersionGenome}/haplotypes.fastq
rm -r ref/${VersionGenome}/primary_index
rm -r ref/${VersionGenome}/tmp
