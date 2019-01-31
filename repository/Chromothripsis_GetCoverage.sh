#!/bin/bash

##########################################################################################
##
## Chromothripsis_GetCoverage.sh
##
## Calculates mean coverage from WGS metrics.
##
##########################################################################################

name=$1

coverage_tumor=$(ps -ef | awk '/MEAN_COVERAGE/ { getline; print $2 }' $name/results/QC/$name.Tumor.bam.metrics)
coverage_normal=$(ps -ef | awk '/MEAN_COVERAGE/ { getline; print $2 }' $name/results/QC/$name.Normal.bam.metrics)
sum=$(echo $coverage_normal + $coverage_tumor | bc)
coverage=$(echo "$sum / 2" | bc)
echo $coverage