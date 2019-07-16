#!/bin/bash

##########################################################################################
##
## config_bash.sh
##
## Adapted for bash.
##
##########################################################################################

repository_dir=/PATH/TO/REPOSITORYDIR
genomes_dir=/PATH/TO/GENOMESDIR

picard_dir=/PATH/TO/PICARDTOOLS
trimmomatic_dir=/PATH/TO/TRIMMOMATIC
GATK_dir=/PATH/TO/GATK
snpeff_dir=/PATH/TO/SNPEFF
hmmcopyutils=/PATH/TO/HMMCOPYUTILS
fasta_to_fastq=/PATH/TO/FASTATOFASTQ
strelka_dir=/PATH/TO/FASTATOSTRELKA
manta_dir=/PATH/TO/FASTATOMANTA
bammatcher_dir=/PATH/TO/BAMMATCHER
vcf2maf_dir=/PATH/TO/VCF2MAF
vep_dir=/PATH/TO/VEP
vepdata_dir=/PATH/TO/VEPDATA

if [ $species = 'Mouse' ]; then
	genome_dir=$genomes_dir/GRCm38.p6
	# only edit fields below if you did not use the initialisation script
	snp_file=$genome_dir/MGP.v5.snp_and_indels.exclude_wild.vcf.gz
	alternate_snp_file=$genome_dir/MGP.v6.snp_and_indels.exclude_wild.vcf.gz
	genome_file=$genome_dir/GRCm38.p6.fna
	genomeindex_dir=$genome_dir/bwa_index/GRCm38.p6
	interval_file=$genome_dir/GRCm38.SureSelect_Mouse_All_Exon_V1.bed.list
	bammatcher_file=$genome_dir/GRCm38.bammatcher_bash.conf
	snpeff_version=GRCm38.86
	microsatellite_file=$genome_dir/GRCm38.p6.microsatellites
	callregions_file=$genome_dir/GRCm38.canonical_chromosomes.bed.gz
	CGC_file=$genome_dir/GRCm38.Census_allMon_Jan_15_11_46_18_2018_mouse.tsv
	TruSight_file="NULL"
	chromosomes=19
	gcWig_file=$genome_dir/GRCm38.p6.gc.20000.wig
	mapWig_file=$genome_dir/GRCm38.p6.map.20000.wig
	exons_file=$genome_dir/GRCm38.SureSelect_Mouse_All_Exon_V1.bed
	centromere_file="NULL"
	varregions_file=$genome_dir/GRCm38.AgilentProbeGaps.txt
	genecode_file=$genome_dir/GRCm38.Genecode_M20_Exons.rds
	vepdata_dir=$genome_dir/VEP
fi