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

cmd_learn_read_orientation="java -jar /opt/gatk-4.1.7.0/gatk.jar LearnReadOrientationModel"
for f in *.${type}.m2.*.f1r2.tar.gz;
do
cmd_learn_read_orientation="${cmd_learn_read_orientation} --input ${f}"
done
cmd_learn_read_orientation="${cmd_learn_read_orientation} --output ${name}.${type}.m2.read-orientation-model.tar.gz"

echo "${cmd_learn_read_orientation}"
${cmd_learn_read_orientation}

