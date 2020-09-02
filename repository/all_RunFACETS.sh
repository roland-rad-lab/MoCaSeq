#!/bin/bash

##########################################################################################
##
## CNV_RunFACETS.sh
##
## Extract SNPs for tumor and normal in order to run FACETS.
##
##########################################################################################

#parallel 'sh '/media/rad/SSD1/DNA/repository/all_RunFACETS.sh' {} Mouse WES /media/rad/SSD1/DNA/configadapted.sh' ::: mPDAC0004_7639_PPT-1

name=$1
species=$2
sequencing_type=$3
config_file=$4

. $config_file

if [ $species = 'Human' ]; then
	#snp_file=$genome_dir/00-common_all.vcf.gz
	snp_file=$dbsnp_file
elif [ $species = 'Mouse' ]; then
	#snp_file=$genome_dir/MGP.v5.snp_and_indels.exclude_wild.chromosomal_sort.vcf.gz
	snp_file=$snp_file # redundant, but like this we know mouse is also covered in this ifelse
fi

if [ $sequencing_type = 'WES' ]; then
	min_coverage="25,0"
elif [ $sequencing_type = 'WGS' ]; then
	min_coverage="15,0"
fi

#rm $name/results/FACETS/$name.pileup.txt.gz

echo -g -q15 -Q20 -P100 -r$min_coverage $snp_file $name/results/FACETS/$name.pileup.txt $name/results/bam/$name.Normal.bam $name/results/bam/$name.Tumor.bam

#~/packages/snp-pileup/./snp-pileup -g -q15 -Q20 -P100 -r$min_coverage $snp_file $name/results/FACETS/$name.pileup.txt $name/results/bam/$name.Normal.bam $name/results/bam/$name.Tumor.bam

/usr/local/lib/R/site-library/facets/extcode/./snp-pileup -g -q15 -Q20 -P100 -r$min_coverage $snp_file $name/results/FACETS/$name.pileup.txt $name/results/bam/$name.Normal.bam $name/results/bam/$name.Tumor.bam
