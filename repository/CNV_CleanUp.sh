#!/bin/bash

##########################################################################################
##
## CNV_CleanUp.sh
##
## Cleans up intermediary files needed for CNV detection.
##
##########################################################################################

name=$1

for chrom in {1..19}; 
do
	rm $name"/results/Copywriter/"$name".SegmentsChromosome"$chrom".txt"
done


