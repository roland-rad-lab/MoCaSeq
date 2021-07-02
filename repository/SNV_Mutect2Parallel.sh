#!/usr/bin/env bash

##########################################################################################
##
## SNV_Mutect2Parallel.sh
##
## Postprocessing for Mutect2.
##
##########################################################################################


name=$1
species=$2
config_file=$3
runmode=$4
artefact_type=$5
GATK=$6
RAM=$7
type=$8

. $config_file

echo "found name '${name}'"
echo "found species '${species}'"
echo "found config_file '${config_file}'"
echo "found runmode '${runmode}'"
echo "found artefact_type '${artefact_type}'"
echo "found GATK '${GATK}'"
echo "found RAM '${RAM}'"
echo "found type '${type}'"

chromosome_names_first=$(echo "${chromosome_names}" | cut -d ' ' -f 1)

mkdir -p $name/results/Mutect2/logs

set -e
# ^ This ensures that the script exits if an internal command fails
# however there are disadvantages, e.g. if grep fails to match anything
# your script will fail. After the parallel section below we turn this off


if [ $runmode = "MS" ]; then
	echo '---- Using Mutect2 Parallel (matched tumor-normal) ----' | tee -a $name/results/QC/$name.report.txt

  for chromosome in ${chromosome_names};
	do
		echo "Run Mutect2 for chromosome ${chromosome}"

		java -Xmx${RAM}G -jar $GATK_dir/gatk.jar Mutect2 \
		--intervals ${chromosome} \
		--native-pair-hmm-threads 2 \
		-R $genome_file \
		-I ${name}/results/bam/${name}.Tumor.bam \
		-I ${name}/results/bam/${name}.Normal.bam \
		-tumor Tumor -normal Normal \
		--f1r2-tar-gz ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.f1r2.tar.gz \
		-O $name/results/Mutect2/${name}.${type}.${chromosome}.m2.vcf \
		-bamout $name/results/Mutect2/${name}.${type}.${chromosome}.m2.bam \
		2> $name/results/Mutect2/logs/${name}.${type}.${chromosome}.log \
		> $name/results/Mutect2/logs/${name}.${type}.${chromosome}.out \
		&
	done
	wait
else
	echo "---- Using Mutect2 Parallel (single-sample) for '${type}' ----" | tee -a $name/results/QC/$name.report.txt

	for chromosome in ${chromosome_names};
	do
		echo "Run Mutect2 for chromosome ${chromosome}"

		java -Xmx${RAM}G -jar $GATK_dir/gatk.jar Mutect2 \
		--intervals ${chromosome} \
		--native-pair-hmm-threads 2 \
		-R $genome_file \
		-I ${name}/results/bam/$name.${type}.bam \
		-tumor $type \
		--f1r2-tar-gz ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.f1r2.tar.gz \
		-O ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.vcf \
		-bamout ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.bam \
		2> ${name}/results/Mutect2/logs/${name}.${type}.${chromosome}.log \
		> ${name}/results/Mutect2/logs/${name}.${type}.${chromosome}.out \
		&
	done
	wait
fi

set +e
# Now we can ignore failing commands again (e.g. grep)

echo '---- Mutect2 Parallel Combine Chromosomes ----' | tee -a ${name}/results/QC/${name}.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

cmd_mutect_stats_merge="java -jar $GATK_dir/gatk.jar MergeMutectStats"
for chromosome in ${chromosome_names};
do
	cmd_mutect_stats_merge="${cmd_mutect_stats_merge} --stats ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.vcf.stats"
done
cmd_mutect_stats_merge="${cmd_mutect_stats_merge} --output ${name}/results/Mutect2/${name}.${type}.m2.vcf.stats"

#echo "${cmd_mutect_stats_merge}"
${cmd_mutect_stats_merge}

# combine VCF files
for chromosome in ${chromosome_names};
do
	if [[ "${chromosome}" == "${chromosome_names_first}" ]]; then
		#echo "chromosome ${chromosome} is first"
		cat ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.vcf > ${name}/results/Mutect2/${name}.${type}.m2.vcf
	else
		#echo "chromosome ${chromosome} rest"
		cat ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.vcf | grep -v "^#" >> ${name}/results/Mutect2/${name}.${type}.m2.vcf
		# grep will return 1 and kill the script if no variants are found || true saves us
	fi
done

if [ $artefact_type = 'yes' ]; then
	cmd_learn_read_orientation="java -jar $GATK_dir/gatk.jar LearnReadOrientationModel"
	for chromosome in ${chromosome_names};
	do
		cmd_learn_read_orientation="${cmd_learn_read_orientation} --input ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.f1r2.tar.gz"
	done
	cmd_learn_read_orientation="${cmd_learn_read_orientation} --output ${name}/results/Mutect2/${name}.${type}.m2.read-orientation-model.tar.gz"

	echo "${cmd_learn_read_orientation}"
	${cmd_learn_read_orientation}
fi

# Cleanup chromosome VCF
for chromosome in ${chromosome_names};
do
	echo "Removing ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.vcf"
	#ls ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.vcf*
	rm ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.vcf*
	rm ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.f1r2.tar.gz
done

# combine BAM files
cmd_samtools_merge="samtools merge -c -p ${name}/results/Mutect2/${name}.${type}.m2.bam"
for chromosome in ${chromosome_names};
do
	cmd_samtools_merge="${cmd_samtools_merge} ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.bam"
done

echo "${cmd_samtools_merge}"
$cmd_samtools_merge

# Cleanup chromosome BAM
for chromosome in ${chromosome_names};
do
	echo "Removing ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.bam"
	rm ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.bam
	rm ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.bai
done

# remove the ".matched" for matched samples (to be consistent with the old format)
if [ $runmode = "MS" ]; then
	for file in ${name}/results/Mutect2/${name}.matched.m2.*; do
			mv "$file" "${file/matched./}"
	done
fi

# Should be safe to assume it's sorted
#samtools sort \
#-@ 4 \
#-o ${name}/results/Mutect2/${name}.${type}.m2.bam \
#-T ${name}/results/Mutect2/${name}.${type}.m2.bam.part. \
#${name}/results/Mutect2/${name}.${type}.m2.unsorted.bam

samtools index ${name}/results/Mutect2/${name}.${type}.m2.bam
#rm ${name}/results/Mutect2/${name}.${type}.m2.unsorted.bam
