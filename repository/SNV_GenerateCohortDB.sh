#!/bin/bash

##########################################################################################
##
## SNV_GenerateCohortDB.sh
##
## Generate a list of all germline positions in > 2 animals from a given cohort.
## Requires (i) list of sample IDs (one per line) and 
## (ii) the raw Mutect2 output (*.NormalOnly.m2.vcf) for all animals in one directory.
##
##########################################################################################

file_list=$1
config_file=$2

. $config_file

parallel --eta --load 80% --noswap -a $file_list 'java -jar '$GATK_dir'/gatk.jar FilterMutectCalls --variant {}.NormalOnly.m2.vcf --output {}.NormalOnly.m2.filt.vcf'

parallel --eta --load 80% --noswap -a $file_list 'cat {}.NormalOnly.m2.filt.vcf | java -jar '$snpeff_dir'/SnpSift.jar filter "( ( na FILTER ) )" > {}.NormalOnly.m2.filt.reduced.vcf'

parallel --eta --load 80% --noswap -a $file_list 'bgzip {}.NormalOnly.m2.filt.reduced.vcf'

parallel --eta --load 80% --noswap -a $file_list 'tabix -p vcf {}.NormalOnly.m2.filt.reduced.vcf.gz'

bcftools isec -c none -n 2+ -O z -o Cohort.txt *NormalOnly.m2.filt.reduced.vcf.gz

printf '##fileformat=VCFv4.2\n' > Cohort.vcf
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n' >> Cohort.vcf
awk -v OFS='\t' '{print $1,$2,".",$3,$4,".","PASS",".","."}' Cohort.txt >> Cohort.vcf