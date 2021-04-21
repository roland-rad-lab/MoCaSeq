#!/bin/bash

##########################################################################################
##
## Chromothripsis_FormatTable.sh
##
## Re-formats rearrangement list provided by Delly.
##
##########################################################################################

name=$1
temp_dir=$2

echo 'Start reformating Tables'

getLines () {
	cut -f 1-2 $1 > part1
	#grep -P -o 'CHR2=[0-9XY]{1,2};POS2=\d*' $1 | sed -r 's/CHR2=([0-9XY]*);END=([0-9]*)/\1\t\2/g' > part2
	grep -P -o 'CHR2=[0-9XY]{1,2}' preciseInsRm.tab | sed -r 's/CHR2=([0-9XY]*)/\1/g' > part2A
	grep -P -o 'POS2=\d*' preciseInsRm.tab | sed -r 's/POS2=([0-9]*)/\1/g' > part2B
	grep -P -o 'PE=\d*;' $1 | sed -r 's/PE=([0-9]*);/\1/g' > part3
	grep -P -o '(IM){0,1}PRECISE' $1 > part5
	cut -f10 $1 | tr : '\t' | cut -f9-12 > part6
	cut -f11 $1 | tr : '\t' | cut -f9-12 > part7
	grep -P -o 'MAPQ=\d*' $1 | sed -r 's/MAPQ=([0-9]*)/\1/g' > part8
	grep -P -o 'CT=\dto\d' $1 | sed -r 's/CT=([35]to[35])/\1/g' > part9
}

mkdir -p ${temp_dir}/${name}_Chromothripsis_tmpData
cd ${temp_dir}/${name}_Chromothripsis_tmpData
inFile=/var/pipeline/$name/results/Delly/$name.pre.vcf

preciseFile='precise.tab'
impreciseFile='imprecise.tab'

#create human-readable translocation table
grep -P '^[0-9XY]{1,2}' $inFile | grep -P 'CHR2=[0-9XY]{1,2}' | grep -P '\bPRECISE' > precise.tab
grep -P '^[0-9XY]{1,2}' $inFile | grep -P 'CHR2=[0-9XY]{1,2}' | grep 'IMPRECISE' > imprecise.tab

#remove Insertion Type varition with NtoN signature
grep -Pv 'INS\d.' precise.tab > preciseInsRm.tab
grep -Pv 'INS\d.' imprecise.tab > impreciseInsRm.tab

#echo 'donChr	donPos	accChr	accPos	PE	SR	Precision	tumor:DR	tumor:DV	tumor:RR	tumor:RV	control:DR	control:DV	control:RR	control:RV	MAPQ	Type' > header
printf 'donChr\tdonPos\taccChr\taccPos\tPE\tSR\tPrecision\ttumor:DR\ttumor:DV\ttumor:RR\ttumor:RV\tcontrol:DR\tcontrol:DV\tcontrol:RR\tcontrol:RV\tMAPQ\tType\n' > header

getLines preciseInsRm.tab
grep -P -o 'SR=\d*;' preciseInsRm.tab | sed -r 's/SR=([0-9]*);/\1/g' > part4
paste part1 part2A part2B part3 part4 part5 part6 part7 part8 part9 > precise2.tab

getLines impreciseInsRm.tab
rm part4; i=0, j=$(wc -l impreciseInsRm.tab);
for i in {1..$j}; do
	echo 'NA' >> part4
done
paste part1 part2A part2B part3 part4 part5 part6 part7 part8 part9 > imprecise2.tab

cat precise2.tab imprecise2.tab > breakpoints2.tab
#cat precise.tab imprecise.tab > breakpoints.tab
sort -k1,1V -k2,2n breakpoints2.tab | cat header - > "/var/pipeline/"$name"/results/Delly/$name".breakpoints.tab

cd /var/pipeline/
rm -r ${temp_dir}/${name}_Chromothripsis_tmpData

echo 'End reformating tables'
