#!/bin/bash

##########################################################################################
##
## CNV_GetGenotype.sh
##
## Extract genotypes for different mouse allels.
##
##########################################################################################

name=$1
position=$2

for type in Normal Tumor;
do
	(
	samtools depth -r $position -a $name/results/bam/$name.$type.bam > $name/results/Genotype/$name.$type.Genotypes.CNV.txt
	) &
done

wait