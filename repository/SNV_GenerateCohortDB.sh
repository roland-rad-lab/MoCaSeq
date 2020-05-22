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

config_file=$1
species=$2

. $config_file

mkdir -p CohortDatabase
mkdir -p CohortDatabase/vcf

parallel --eta --load 80% --noswap \
 'java -jar '$GATK_dir'/gatk.jar FilterMutectCalls \
 --variant {}/results/Mutect2/{}.Normal.m2.vcf \
 --output CohortDatabase/vcf/{}.Normal.m2.filt.vcf \
 --reference '$genome_file'' ::: $(ls -l | grep "^d" | grep "PPT" | awk '{print $9}')

parallel --eta --load 80% --noswap \
'cat CohortDatabase/vcf/{}.Normal.m2.filt.vcf \
| java -jar '$snpeff_dir'/SnpSift.jar filter \
"( ( FILTER = '\''PASS'\'' ) | ( FILTER = '\''clustered_events'\'' ) | ( FILTER = '\''str_contraction'\'' ) )" \
> CohortDatabase/vcf/{}.Normal.m2.filt.reduced.vcf; \
rm CohortDatabase/vcf/{}.Normal.m2.filt.vcf; \
rm CohortDatabase/vcf/{}.Normal.m2.filt.vcf.idx' ::: $(ls -l | grep "^d" | grep "PPT" | awk '{print $9}')

parallel --eta --load 80% --noswap \
'cat CohortDatabase/vcf/{}.Normal.m2.filt.reduced.vcf \
| java -jar '$snpeff_dir'/SnpSift.jar filter \
"( (GEN[Normal].AD[1] >= 2) )" \
> CohortDatabase/vcf/{}.Normal.m2.filt.reduced.read_filter.vcf; \
rm CohortDatabase/vcf/{}.Normal.m2.filt.reduced.vcf; \
rm CohortDatabase/vcf/{}.Normal.m2.filt.vcf.reduced.idx' ::: $(ls -l | grep "^d" | grep "PPT" | awk '{print $9}')

parallel --eta --load 80% --noswap \
'bgzip CohortDatabase/vcf/{}.Normal.m2.filt.reduced.read_filter.vcf' ::: $(ls -l | grep "^d" | grep "PPT" | awk '{print $9}')

parallel --eta --load 80% --noswap \
'tabix -p vcf CohortDatabase/vcf/{}.Normal.m2.filt.reduced.read_filter.vcf.gz' ::: $(ls -l | grep "^d" | grep "PPT" | awk '{print $9}')

bcftools isec -c none -n 2+ -O z -o CohortDatabase/Cohort.txt CohortDatabase/vcf/*.Normal.m2.filt.reduced.read_filter.vcf.gz

printf '##fileformat=VCFv4.2\n' > CohortDatabase/Cohort.vcf
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n' >> CohortDatabase/Cohort.vcf
awk -v OFS='\t' '{print $1,$2,".",$3,$4,".","PASS",".","."}' CohortDatabase/Cohort.txt >> CohortDatabase/Cohort.vcf

bgzip CohortDatabase/Cohort.vcf
tabix -p vcf CohortDatabase/Cohort.vcf.gz

bcftools stats CohortDatabase/Cohort.vcf.gz > CohortDatabase/Cohort.vcf.gz.stats

bcftools norm -m -any CohortDatabase/Cohort.vcf.gz \
-O z -o CohortDatabase/Cohort.fixed.vcf.gz

tabix -p vcf CohortDatabase/Cohort.fixed.vcf.gz

bcftools stats CohortDatabase/Cohort.fixed.vcf.gz > CohortDatabase/Cohort.fixed.vcf.gz.stats

mv CohortDatabase/Cohort.fixed.vcf.gz CohortDatabase/Cohort.vcf.gz
mv CohortDatabase/Cohort.fixed.vcf.gz.tbi CohortDatabase/Cohort.vcf.gz.tbi

rm CohortDatabase/Cohort.txt