#!/bin/bash

##########################################################################################
##
## CNV_GetGenotype.sh
##
## Extract genotypes for different mouse allels.
##
##########################################################################################

name=$1

for chrom in $(seq 22);
do
	rm $name"/results/Copywriter/"$name".SegmentsChromosome"$chrom".txt"
done