#!/bin/bash
type=-1
species=Human
config_file=/opt/MoCaSeq/config.sh
filtering="soft"
artefact_type="yes"
GATK=4.1.7.0
repository_dir=/opt/MoCaSeq/repository

# parse parameters
if [ "$1" = "" ]; then usage; fi
while [ "$1" != "" ]; do case $1 in
	-t|--type) shift;type="$1";;
    --help) usage;shift;;
	*) usage;shift;;
esac; shift; done

# check if parameter type was specified
if [[ $type == -1 ]]; then echo "no type specified!" &&	exit 1; fi

. $config_file

if [[ $type == "matched" ]]
then
	echo "detected type $type"

	bash $repository_dir/SNV_Mutect2Postprocessing_LRZ.sh $name $species $config_file $filtering $artefact_type $GATK $type

	Rscript $repository_dir/SNV_SelectOutput.R $name Mutect2 $species $CGC_file $TruSight_file
else
	echo "detected type $type"
	
	bash $repository_dir/SNV_Mutect2PostprocessingSS_LRZ.sh $name $species $config_file $filtering $artefact_type $GATK $type

	Rscript $repository_dir/SNV_SelectOutputSS.R $name Mutect2 $species $CGC_file $TruSight_file
fi

Rscript $repository_dir/SNV_Signatures.R $name $species

