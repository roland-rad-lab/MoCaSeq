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
  file_to_encrypt=$1

    if [ -n "$(find "$file_to_encrypt.c4gh" -prune -size +1c)" ]; then
        echo "File successfully encrypted: ${file_to_encrypt}"
        rm ${file_to_encrypt}
    else
        echo "File was NOT successfully encrypted: ${file_to_encrypt}"
        echo "This file was NOT encrypted! (maybe the password in argument 2 is wrong?)"
        rm ${file_to_encrypt}.c4gh # remove empty file 
    fi
}

function encrypt{
  file_to_encrypt=$1
  
  if [ -f $file_to_encrypt ]
  then
  encrypt $file_to_encrypt
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
## Normal bam file
file_to_encrypt=${name}/results/bam/${name}.Normal.bam
encrypt $file_to_encrypt

## Normal Mutect files
file_to_encrypt=${name}/results/Mutect2/${name}.Normal.m2.bam
encrypt $file_to_encrypt

file_to_encrypt=${name}/results/Mutect2/${name}.Normal.m2.vcf
encrypt $file_to_encrypt

file_to_encrypt=${name}/results/Mutect2/${name}.Normal.Mutect2.txt
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Normal.Mutect2.vcf
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Normal.Mutect2.vep.maf
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Normal.Mutect2.vep.vcf
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Normal.Mutect2.mergeid.vcf
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Normal.m2.filt.selected.vcf
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Normal.Mutect2.Positions.txt
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Normal.Mutect2.NoCommonSNPs.txt
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Normal.Mutect2.annotated.one.vcf
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Normal.Mutect2.NoCommonSNPs.OnlyImpact.txt
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Normal.Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Normal.Mutect2.NoCommonSNPs.OnlyImpact.TruSight.txt
encrypt $file_to_encrypt


elif [[ $type == "Tumor" ]]

# TUMOR file encryption
## Tumor bam file
file_to_encrypt=${name}/results/bam/${name}.Tumor.bam
encrypt $file_to_encrypt


# LOH
file_to_encrypt=${name}/results/LOH/${name}.VariantsForLOH.txt
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/LOH/${name}.VariantsForLOHGermline.txt
encrypt $file_to_encrypt


# PURITY
file_to_encrypt=${name}/results/InhousePurity/${name}_corrected_LOH_values.tsv.gz
encrypt $file_to_encrypt

file_to_encrypt=${name}/results/InhousePurity/Flex/${name}_corrected_LOH_values.tsv.gz
encrypt $file_to_encrypt

# ALL THE MUTECT2 FILES

# COMBINED
file_to_encrypt=${name}/results/Mutect2/${name}.m2.bam
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.m2.vcf
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.m2.filt.vcf
encrypt $file_to_encrypt


# TUMOR
file_to_encrypt=${name}/results/Mutect2/${name}.Tumor.m2.bam
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Tumor.m2.vcf
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Tumor.Mutect2.txt
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Tumor.Mutect2.vcf
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Tumor.Mutect2.vep.maf
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Tumor.Mutect2.vep.vcf
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Tumor.Mutect2.mergeid.vcf
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Tumor.m2.filt.selected.vcf
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Tumor.Mutect2.Positions.txt
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Tumor.Mutect2.NoCommonSNPs.txt
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Tumor.Mutect2.annotated.one.vcf
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Tumor.Mutect2.NoCommonSNPs.OnlyImpact.txt
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Tumor.Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt
encrypt $file_to_encrypt


file_to_encrypt=${name}/results/Mutect2/${name}.Tumor.Mutect2.NoCommonSNPs.OnlyImpact.TruSight.txt
encrypt $file_to_encrypt


else
echo "invalid type: $type"
exit 1

fi

