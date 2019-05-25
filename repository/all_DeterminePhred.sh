#!/bin/bash

##########################################################################################
##
## all_DeterminePhred.sh
##
## Extracts Phred-Score from FASTQ-output.
##
##########################################################################################

name=$1
type=$2

phred=$(unzip -p $name/results/QC/$name.$type.R1_fastqc.zip $name.$type.R1_fastqc/fastqc_data.txt | grep Encoding | awk '{print $5}')

if [ "$phred" = "1.9" ] || [ "$phred" = "1.8" ]; then
	phred=phred33
else phred=phred64
fi

echo $phred