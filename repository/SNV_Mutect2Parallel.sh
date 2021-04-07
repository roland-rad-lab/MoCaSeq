#!/bin/bash

##########################################################################################
##
## SNV_Mutect2Parallel.sh
##
## Postprocessing for Mutect2.
##
##########################################################################################

set -e
# ^ This ensures that the script exits if an internal command fails

chromosome_names="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y"
chromosome_names_first=$(echo "${chromosome_names}" | cut -d ' ' -f 1)

name=$1
species=$2
config_file=$3
runmode=$4
artefact_type=$5
GATK=$6
type=$7
# ^No type unless runmode == SS (single sample)

. $config_file

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
		2> $name/results/Mutect2/$name.${type}.${chromosome}.log \
		> $name/results/Mutect2/$name.${type}.${chromosome}.out \
		&
	done
	wait
else
	echo '---- Using Mutect2 Parallel (single-sample) ----' | tee -a $name/results/QC/$name.report.txt

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
		2> ${name}/results/Mutect2/${name}.${type}.${chromosome}.log \
		> ${name}/results/Mutect2/${name}.${type}.${chromosome}.out \
		&
	done
	wait
fi

echo '---- Mutect2 Parallel Combine Chromosomes ----' | tee -a ${name}/results/QC/${name}.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

cmd_mutect_stats_merge="java -jar $GATK_dir/gatk.jar MergeMutectStats"
for chromosome in ${chromosome_names};
do
	cmd_mutect_stats_merge="${cmd_mutect_stats_merge} --stats ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.vcf.stats"
done
cmd_mutect_stats_merge="${cmd_mutect_stats_merge} --output ${name}/results/Mutect2/${name}.${type}.m2.combined.vcf.stats"

#echo "${cmd_mutect_stats_merge}"
${cmd_mutect_stats_merge}

# combine VCF files
for chromosome in ${chromosome_names};
do
	if [[ $chromosome == $chromosome_names_first ]]; then
		#echo "chromosome ${chromosome} is first"
		cat ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.vcf > ${name}/results/Mutect2/${name}.${type}.m2.vcf
	else
		#echo "chromosome ${chromosome} rest"
		cat ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.vcf | grep -v "^#" >> ${name}/results/Mutect2/${name}.${type}.m2.vcf
	fi
done

# combine BAM files
cmd_samtools_merge="samtools merge ${name}/results/Mutect2/${name}.${type}.m2.unsorted.bam"
for chromosome in ${chromosome_names};
do
	cmd_samtools_merge="${cmd_samtools_merge} ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.bam"
done

echo "${cmd_samtools_merge}"
$cmd_samtools_merge

for chromosome in ${chromosome_names};
do
	echo "Removing ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.bam"
	rm ${name}/results/Mutect2/${name}.${type}.${chromosome}.m2.bam
done

samtools sort \
-@ 4 \
-o ${name}/results/Mutect2/${name}.${type}.m2.bam \
-T ${name}/results/Mutect2/${name}.${type}.m2.bam.part. \
${name}/results/Mutect2/${name}.${type}.m2.unsorted.bam

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
