#!/bin/bash

##########################################################################################
##
## LOH_CNVKitPrepareLOH.sh
##
## Prepare LOH positions for CNVKit.
##
##########################################################################################

name=$1
genome_file=$2

cat $name/results/LOH/$name.VariantsForLOH.txt | cut -f2,3 | tail -n +2 > $name/results/LOH/$name.VariantsForLOH.sites.txt

bcftools mpileup --ignore-RG --skip-indels --count-orphans -a AD,DP \
-f $genome_file -R $name/results/LOH/$name.VariantsForLOH.sites.txt \
$name/results/bam/$name.Tumor.bam | bcftools call - -O v -mv -A \
-o $name/results/LOH/$name.VariantsForLOH.sites.vcf