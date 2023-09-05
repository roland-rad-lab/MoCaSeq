#!/bin/bash

name=$1
password=$2
type=$3 # Tumor or Normal

echo "Sample: $name";
export C4GH_PASSPHRASE=$password

# define the path to the pass files
skFile=~/c4ghKeys/COMPASS_c4gh.sec
pkFile=~/c4ghKeys/COMPASS_c4gh.pub

# function to check if a file was successfully encrypted or not
# check if c4gh is valid (i.e. larger than 1 Byte) (i.e. not delete if the given password was wrong)
function CheckEncryption {
  file_to_check=$1

    if [ -n "$(find "$file_to_check.c4gh" -prune -size +1c)" ]; then
        echo "File successfully encrypted: ${file_to_check}"
        rm ${file_to_check}
    else
        echo "File was NOT successfully encrypted: ${file_to_check}"
        echo "This file was NOT encrypted! (maybe the password in argument 2 is wrong?)"
        rm ${file_to_check}.c4gh # remove empty file 
    fi
}

function encrypt {
  file_to_encrypt=$1
  
  if [ -f $file_to_encrypt ]
  then
  ~/.local/bin/crypt4gh encrypt --sk $skFile --recipient_pk $pkFile < ${file_to_encrypt} > ${file_to_encrypt}.c4gh
  CheckEncryption $file_to_encrypt
  else
  echo "File NOT found: $file_to_encrypt"
  fi
}


# ENCRYPT MULTIPLE SPECIFIC FILES 

# file_to_encrypt=XXXXXXXXX
# encrypt $file_to_encrypt

if [[ $type == "Normal" ]]
then

# NORMAL file encryption
# bam files
file_to_encrypt=${name}/results/bam/${name}.Normal.bam
encrypt $file_to_encrypt
file_to_encrypt=${name}/results/Mutect2/${name}.Normal.m2.bam
encrypt $file_to_encrypt


# vcfs
# vcfs
vcf_files=(${name}/results/Mutect2/${name}.Normal.m2.vcf.gz
${name}/results/Mutect2/${name}.Normal.m2.filtered.vcf.gz 
${name}/results/Mutect2/${name}.Normal.Mutect2.vcf.gz
${name}/results/Mutect2/${name}.Normal.Mutect2.vep.vcf
${name}/results/Mutect2/${name}.Normal.Mutect2.annotated.vcf.gz)
# these files are not created by the nextflow version of MoCaSeq
# ${name}/results/Mutect2/${name}.Normal.Mutect2.mergeid.vcf
# ${name}/results/Mutect2/${name}.Normal.m2.filt.selected.vcf
# ${name}/results/Mutect2/${name}.Normal.Mutect2.annotated.one.vcf -> ${name}/results/Mutect2/${name}.*.Mutect2.annotated.vcf.gz

for file_to_encrypt in ${vcf_files[@]};
do
encrypt $file_to_encrypt;
done

# custom vcf style txt files
txt_files=(${name}/results/Mutect2/${name}.Normal.Mutect2.txt
${name}/results/Mutect2/${name}.Normal.Mutect2.Positions.txt
${name}/results/Mutect2/${name}.Normal.Mutect2.NoCommonSNPs.txt
${name}/results/Mutect2/${name}.Normal.Mutect2.NoCommonSNPs.OnlyImpact.txt
${name}/results/Mutect2/${name}.Normal.Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt
${name}/results/Mutect2/${name}.Normal.Mutect2.NoCommonSNPs.OnlyImpact.TruSight.txt)

for file_to_encrypt in ${txt_files[@]};
do
encrypt $file_to_encrypt;
done


elif [[ $type == "Tumor" ]]
then
# TUMOR file encryption
## bam files
bam_files=(${name}/results/bam/${name}.Tumor.bam 
${name}/results/Mutect2/${name}.m2.bam
${name}/results/Mutect2/${name}.Tumor.m2.bam)

for file_to_encrypt in ${bam_files[@]};
do
encrypt $file_to_encrypt;
done


# LOH
file_to_encrypt=${name}/results/LOH/${name}.VariantsForLOH.txt
encrypt $file_to_encrypt

file_to_encrypt=${name}/results/LOH/${name}.VariantsForLOHGermline.txt
encrypt $file_to_encrypt


# PURITY
# these are pre-deleted by Niklas
#${name}/results/InhousePurity/${name}_corrected_LOH_values.tsv.gz
#${name}/results/InhousePurity/Flex/${name}_corrected_LOH_values.tsv.gz

# vcfs
vcf_files=(${name}/results/Mutect2/${name}.*.m2.vcf.gz
${name}/results/Mutect2/${name}.*.m2.filtered.vcf.gz 
${name}/results/Mutect2/${name}.*.Mutect2.vcf.gz
${name}/results/Mutect2/${name}.*.Mutect2.vep.vcf
${name}/results/Mutect2/${name}.*.Mutect2.annotated.vcf.gz)

for file_to_encrypt in ${vcf_files[@]};
do
encrypt $file_to_encrypt;
done

# custom vcf style txt files
txt_files=(${name}/results/Mutect2/${name}.*.Mutect2.txt
${name}/results/Mutect2/${name}.*.Mutect2.Positions.txt
${name}/results/Mutect2/${name}.Tumor.Mutect2.NoCommonSNPs.txt
${name}/results/Mutect2/${name}.Tumor.Mutect2.NoCommonSNPs.OnlyImpact.txt
${name}/results/Mutect2/${name}.Tumor.Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt
${name}/results/Mutect2/${name}.Tumor.Mutect2.NoCommonSNPs.OnlyImpact.TruSight.txt)

for file_to_encrypt in ${txt_files[@]};
do
encrypt $file_to_encrypt;
done

# maf
file_to_encrypt=${name}/results/Mutect2/${name}.Tumor.Mutect2.vep.maf
encrypt $file_to_encrypt

else
echo "Invalid type: $type"
exit 1

fi

