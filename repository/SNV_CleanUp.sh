#!/bin/bash

##########################################################################################
##
## SNV_CleanUp.sh
##
## Removes all intermediary files from SNV/Indel and LOH detection.
##
##########################################################################################

name=$1
analysis=$2

if [ $analysis = 'MS' ]; then
	rm $name/results/Mutect2/$name.m2.filt.AM.filtered.selected.vcf.idx
	rm $name/results/Mutect2/$name.m2.filt.AM.filtered.selected.vcf
	rm $name/results/Mutect2/$name.m2.filt.AM.filtered.vcf
	rm $name/results/Mutect2/$name.m2.filt.AM.filtered.vcf.idx
	rm $name/results/Mutect2/$name.m2.filt.AM.vcf
	rm $name/results/Mutect2/$name.m2.filt.AM.vcf.idx
	rm $name/results/Mutect2/$name.m2.filt.AM.vcf.summary
	rm $name/results/Mutect2/$name.m2.filt.vcf
	rm $name/results/Mutect2/$name.m2.filt.vcf.idx
	rm $name/results/Mutect2/$name.m2.postprocessed.vcf.gz
	rm $name/results/Mutect2/$name.m2.postprocessed.vcf.gz.tbi
	rm $name/results/Mutect2/$name.Mutect2.annotated.one.vcf
	rm $name/results/Mutect2/$name.Mutect2.annotated.vcf
	rm $name/results/Mutect2/$name.Mutect2.annotated.vcf.stats
	rm $name/results/Mutect2/$name.Mutect2.annotated.vcf.stats.genes.txt
	rm snpEff_summary.html
elif [ $analysis = 'SS' ]; then
	for type in Normal Tumor;
	do
	rm $name/results/Mutect2/$name.$type.m2.filt.AM.reduced.vcf
	rm $name/results/Mutect2/$name.$type.m2.filt.AM.vcf
	rm $name/results/Mutect2/$name.$type.m2.filt.AM.vcf.idx
	rm $name/results/Mutect2/$name.$type.m2.filt.AM.vcf.summary
	rm $name/results/Mutect2/$name.$type.m2.filt.vcf
	rm $name/results/Mutect2/$name.$type.m2.filt.vcf.idx
	rm $name/results/Mutect2/$name.$type.m2.vcf.idx
	done
fi