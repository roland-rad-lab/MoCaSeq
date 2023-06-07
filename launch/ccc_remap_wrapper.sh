#!/bin/bash
# charliecloud remapping wrapper MoCaSeq pipeline
# @author:marcus.wagner@tum.de
# 

usage()
{
	echo "  Usage: $0 "
	echo "	-ccc	--charliecloud_continer	Path to charliecloud continer to use."
	echo "	-wd	--working_directory	Path to working directory."
	echo "	-m	--mocaseq		Path to MoCaSeq git repository clone."
	echo "	-s	--sample		Sample name, used for naming output and intermediate files."
	echo "	-bd	--bam_dir		Path to bam input file for remapping."
	echo "	-bf	--bam_file		.bam input file name."
	echo "	-t	--tumor			Tumor sample flag."
	echo "	-h	--help			Display this message."
  exit 1
}

# default parameters
cccDir=
workingDir=
mocaseqDir=
bamDir=
bamName=
sample=
type="normal"

# parse parameters
if [ "$1" = "" ]; then usage; fi
while [ "$1" != "" ]; do case $1 in
	-ccc|--charliecloud_continer) shift;cccDir="$1";;
	-s|--sample) shift;sample="$1";;
	-bf|--bam_file) shift;bamName="$1";;
	-bd|--bam_dir) shift;bamDir="$1";;
	-wd|--working_directory) shift;workingDir="$1";;
	-m|--mocaseq) shift;mocaseqDir="$1";;
	-t|--tumor) shift;type="tumor";;
    --help) usage;shift;;
	*) usage;shift;;
esac; shift; done

# charliecloud call
if [[ $type == "normal" ]]
then
	ch-run $cccDir --no-home --set-env=sample=${sample} -w --no-passwd \
--bind ${workingDir}:/var/pipeline/ \
--bind ${mocaseqDir}:/opt/MoCaSeq/ \
--bind ${bamDir}:/var/raw-bams/ \
-- \
/opt/MoCaSeq/MoCaSeq_LRZ_remap.sh \
-nb /var/raw-bams/${bamName} \
--name ${sample} \
--species Human \
--repeat_mapping yes \
--sequencing_type WGS \
--quality_control yes \
--threads 30 \
--RAM 45 \
--GATKVersion 4.1.7.0 \
--filtering soft \
--artefact yes
else
	ch-run $cccDir --no-home --set-env=sample=${sample} -w --no-passwd \
--bind ${workingDir}:/var/pipeline/ \
--bind ${mocaseqDir}:/opt/MoCaSeq/ \
--bind ${bamDir}:/var/raw-bams/ \
-- \
/opt/MoCaSeq/MoCaSeq_LRZ_remap.sh \
-tb /var/raw-bams/${bamName} \
--name ${sample} \
--species Human \
--repeat_mapping yes \
--sequencing_type WGS \
--quality_control yes \
--threads 30 \
--RAM 45 \
--GATKVersion 4.1.7.0 \
--filtering soft \
--artefact yes
fi

