#!/bin/bash
name=$1
species=Human
config_file=/opt/MoCaSeq/config.sh
filtering="soft"
artefact_type="yes"
GATK=4.1.7.0
type=matched
repository_dir=/opt/MoCaSeq/repository

. $config_file

bash $repository_dir/SNV_Mutect2Postprocessing_LRZ.sh $name $species $config_file \
 $filtering $artefact_type $GATK $type

Rscript $repository_dir/SNV_SelectOutput.R $name Mutect2 $species $CGC_file $TruSight_file

Rscript $repository_dir/SNV_Signatures.R $name $species

