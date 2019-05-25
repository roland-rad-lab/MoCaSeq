#!/bin/bash

##########################################################################################
##
## SNV_GetGenotype.sh
##
## Extract genotypes for different mouse allels.
##
##########################################################################################

name=$1
allele=$2
comment=$3
config_file=$4
species=$5
position=$6
runmode=$7
types=$8

. $config_file

chrom=$(echo $position | cut -f1 -d':')
start=$(echo $position | cut -f1 -d'-' | cut -f2 -d':')
end=$(echo $position | cut -f2 -d'-')

if [ $runmode = "MS" ]; then
	types="Tumor Normal"
fi

for type in $types;
do
	samtools mpileup -A -B -f $genome_file -q 0 -Q 20 -v -u -t AD \
	$name/results/bam/$name.$type.bam -r $position | \
	java -jar $snpeff_dir/SnpSift.jar extractFields - \
	CHROM POS REF ALT FORMAT[0].AD[0] FORMAT[0].AD[1] \
	 > $name/results/Genotype/$name.Genotypes.temp.SNV.txt

	var=$(wc -l < $name/results/Genotype/$name.Genotypes.temp.SNV.txt)

	if [ "$var" -gt "1" ]; then
		awk 'BEGIN {OFS="\t"}; FNR>1 {print "'$name'_'$type'","'$allele'",$0,"'$comment'"}' \
		 $name/results/Genotype/$name.Genotypes.temp.SNV.txt  \
		 >> $name/results/Genotype/$name.Genotypes.txt
		 rm $name/results/Genotype/$name.Genotypes.temp.SNV.txt
	elif [ "$var" -eq "1"  ]; then
		echo ''$name'_'$type'\t'$allele'\t'$chrom'\t'$start'\t'$end'\tNoReads\t\t\t' \
		 >> $name/results/Genotype/$name.Genotypes.txt
	fi
done