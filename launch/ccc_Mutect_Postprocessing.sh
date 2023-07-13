#!/bin/bash

bash $repository_dir/SNV_Mutect2Postprocessing_LRZ.sh $name $species $config_file \
$name $species $config_file $filtering $artefact_type $GATK $type

Rscript $repository_dir/SNV_SelectOutput.R $name Mutect2 $species $CGC_file $TruSight_file

Rscript $repository_dir/SNV_Signatures.R $name $species

