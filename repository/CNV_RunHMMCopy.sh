#!/bin/bash

##########################################################################################
##
## CNV_RunHMMCopy.sh
##
## Run HMMCopy und matched tumour-normal .bam-files.
##
##########################################################################################

name=$1
species=$2
config_file=$3
resolution=$4 # window size for the CNV estimation; usually between 1000 and 20000 (depending on the sequencing coverage)
types=$5
additionalType=$6 # just used to catch invalid $type parameter (no quotes)

# in the publication the variables are given as: sh $name Mouse $config_file 20000 Tumor Normal
# the separated naming of Tumor and Normal (no quotes) is an unsupported format, but this dirty fix will help
# check and fix if there is a "unexpected" value given after types --> which results from using Tumor Normal or $types instead of "Tumor Normal"
if [ "$additionalType" = "Normal" ]; then
	types="Tumor Normal"
fi

. $config_file

chromosomes=$(echo "${chromosome_names}" | tr " " ",")

for type in $types;
do
	(
	echo "Binning read counts in $type file @ $resolution resolution..."
	echo "Binning read counts for ${chromosomes}..."
	# readCounter requires indices as .bam.bai (there is a -b option to build them automatically but it doesn't check if they exist first)
	$hmmcopyutils_dir/bin/readCounter -w $resolution -q20 -c $chromosomes $name/results/bam/$name.$type.bam > $name/results/HMMCopy/$name.$type.$resolution.wig
	) &
done

wait
