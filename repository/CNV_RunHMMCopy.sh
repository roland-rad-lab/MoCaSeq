#!/bin/bash

##########################################################################################
##
## CNV_RunHMMCopy.sh
##
## Run HMMCopy und matched tumour-normal .bam-files.
##
##########################################################################################

name=$1
resolution=$2 # window size for the CNV estimation; usually between 1000 and 20000 (depending on the sequencing coverage)

echo "Binning read counts in control file @ $resolution resolution..."
~/packages/hmmcopy_utils-1.0/bin/readCounter -w $resolution -q20 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y $name/results/bam/$name.Normal.bam > $name/results/HMMCopy/$name.Normal.$resolution.wig
echo "Binning read counts in tumor file @ $resolution resolution..."
~/packages/hmmcopy_utils-1.0/bin/readCounter -w $resolution -q20 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y $name/results/bam/$name.Tumor.bam > $name/results/HMMCopy/$name.Tumor.$resolution.wig