#!/bin/bash

##########################################################################################
##
## all_MoCaSeq_lcWGS.sh
##
## Main workflow for lcWGS
##
##########################################################################################

usage()
{
	echo "  Usage: $0 "
	echo "	-n, --name               Name of the sample."
	echo "	-nf, --fastq_normal_fw   Path to first normal Fastq. Do NOT use if running single-sample tumor only."
	echo "	-nr, --fastq_normal_rev  Path to second normal Fastq. Do NOT use if running single-sample tumor only."
	echo "	-tf, --fastq_tumor_fw    Path to first tumor fastq. Do NOT use if running single-sample normal only."
	echo "	-tr, --fastq_tumor_rev   Path to second tumor fastq. Do NOT use if running single-sample normal only."
	echo "	-e, --ends               Determine sequencing mode. Choose from PE or SE. Defaults to SE."	
	echo "	-t, --threads            Number of CPU threads. Optional. Defaults to 8."
	echo "	-r, --RAM                Amount of Gb RAM. Optional. Defaults to 32."
	echo "	-temp, --temp_dir        Path to temporary directory. Optional. Defaults to current working directory."
	echo "	-p, --phred              If not set, script will try to automatically extract phred-score. Otherwise, set manually to 'phred33' or 'phred64'. 'phred64' only relevant for Illumina data originating before 2011. Optional."
	echo "	--help                   Show this help."
  exit 1
}

# default parameters
fastq_normal_1=
fastq_normal_2=
fastq_tumor_1=
fastq_tumor_2=
sequencing_type=lcWGS
species=Mouse
quality_control=yes
threads=8
RAM=32
temp_dir=/var/pipeline/temp
phred=phred33
runmode=MS
types="Tumor Normal"
config_file=
GATK=4.1.3.0
chromosomes=19
ends=SE

# parse parameters
if [ "$1" = "" ]; then usage; fi
while [ "$1" != "" ]; do case $1 in
	-n|--name) shift;name="$1";;
	-nf|--fastq_normal_fw) shift;fastq_normal_1="$1";;
	-nr|--fastq_normal_rev) shift;fastq_normal_2="$1";;
	-tf|--fastq_tumor_fw) shift;fastq_tumor_1="$1";;
	-tr|--fastq_tumor_rev) shift;fastq_tumor_2="$1";;
	-e|--ends) shift;ends="$1";;
	-t|--threads) shift;threads="$1";;
	-r|--RAM) shift;RAM="$1";;
	-temp|--temp_dir) shift;temp_dir="$1";;
    -p|--phred) shift;phred="$1";;
    --test) shift;test="$1";;
    --help) usage;shift;;
	*) usage;shift;;
esac; shift; done

config_file=/opt/MoCaSeq/config.sh

#reading configuration from $config_file
source $config_file
repository_dir=${config_file%/*}/repository

echo '---- Starting Mouse Cancer Genome Analysis ----'
echo -e "$(date) \t timestamp: $(date +%s)"

echo '---- Creating directories ----'
echo -e "$(date) \t timestamp: $(date +%s)"
mkdir -p $name/
mkdir -p $name/pipeline/
mkdir -p $name/fastq/
mkdir -p $name/results/QC
mkdir -p $name/results/bam
mkdir -p $name/results/Copywriter
mkdir -p $name/results/HMMCopy

if [ ! -d $temp_dir ]; then
  mkdir -p $temp_dir/
fi

if [ $RAM -ge 16 ]; then
	bwainputbases=100000000
else bwainputbases=10000000
fi

MAX_RECORDS_IN_RAM=$(expr $RAM \* 250000)
HASH_TABLE_SIZE=$((RAM*1000000000/500))

echo '---- Starting Mouse Cancer Genome Analysis ----' | tee -a $name/results/QC/$name.report.txt
echo Starting lcWGS pipeline using these settings: | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
echo Running sample named $name | tee -a $name/results/QC/$name.report.txt
echo Running in $ends-mode | tee -a $name/results/QC/$name.report.txt
if [ $runmode = "MS" ] && [ ! -z $fastq_normal_1 ]; then
	echo Using $fastq_normal_1 and $fastq_normal_2 for normal fastqs | tee -a $name/results/QC/$name.report.txt
	echo Using $fastq_tumor_1 and $fastq_tumor_2 for tumor fastqs | tee -a $name/results/QC/$name.report.txt
elif [ $runmode = "SS" ] && [ $repeat_mapping = "no" ]; then
	echo Assigning $fastq_tumor_1 and $fastq_tumor_2 as $types | tee -a $name/results/QC/$name.report.txt
elif [ $runmode = "SS" ] && [ $repeat_mapping = "no" ]; then
	echo Assigning $fastq_normal_1 and $fastq_normal_2 as $types | tee -a $name/results/QC/$name.report.txt
fi

echo Assuming that reads are from $species | tee -a $name/results/QC/$name.report.txt
echo Assuming that experiment is $sequencing_type | tee -a $name/results/QC/$name.report.txt
echo Reading configuration file from $config_file | tee -a $name/results/QC/$name.report.txt
echo Setting location of repository to $repository_dir | tee -a $name/results/QC/$name.report.txt
echo Setting location of genome to $genome_dir | tee -a $name/results/QC/$name.report.txt
echo Setting location for temporary files to $temp_dir| tee -a $name/results/QC/$name.report.txt
echo Quality scores are assumed as $phred | tee -a $name/results/QC/$name.report.txt
echo Starting workflow using $threads CPU-threads and $RAM GB of RAM | tee -a $name/results/QC/$name.report.txt

#rerouting STDERR to report file
exec 2>> $name/results/QC/$name.report.txt

echo '---- Creating directories ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

echo '---- Copying repository ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

cp -r $repository_dir/ $name/pipeline/
cp $repository_dir/../MoCaSeq.sh $name/pipeline/
cp $config_file $name/pipeline/

echo '---- Copying raw data ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

if [ $runmode = 'MS' ] && [ $ends = 'PE' ]; then
cp $fastq_normal_1 $name/fastq/$name.Normal.R1.fastq.gz
cp $fastq_normal_2 $name/fastq/$name.Normal.R2.fastq.gz
cp $fastq_tumor_1 $name/fastq/$name.Tumor.R1.fastq.gz
cp $fastq_tumor_2 $name/fastq/$name.Tumor.R2.fastq.gz

elif [ $runmode = 'SS' ] && [ $types = 'Tumor' ] && [ $ends = 'PE' ]; then
cp $fastq_tumor_1 $name/fastq/$name.$types.R1.fastq.gz
cp $fastq_tumor_2 $name/fastq/$name.$types.R2.fastq.gz

elif [ $runmode = 'SS' ] && [ $types = 'Normal' ] && [ $ends = 'PE' ]; then
cp $fastq_normal_1 $name/fastq/$name.$types.R1.fastq.gz
cp $fastq_normal_2 $name/fastq/$name.$types.R2.fastq.gz

elif [ $runmode = 'MS' ] && [ $ends = 'SE' ]; then
cp $fastq_normal_1 $name/fastq/$name.Normal.R1.fastq.gz
cp $fastq_tumor_1 $name/fastq/$name.Tumor.R1.fastq.gz

elif [ $runmode = 'SS' ] && [ $types = 'Tumor' ] && [ $ends = 'SE' ]; then
cp $fastq_tumor_1 $name/fastq/$name.$types.R1.fastq.gz

elif [ $runmode = 'SS' ] && [ $types = 'Normal' ] && [ $ends = 'SE' ]; then
cp $fastq_normal_1 $name/fastq/$name.$types.R1.fastq.gz
fi

if [ $ends = 'PE' ]; then
	echo '---- Calculating md5-sums ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
	for type in $types;
	do
		md5sum $name/fastq/$name.$type.R1.fastq.gz > $name/fastq/$name.$type.R1.fastq.gz.md5
		md5sum $name/fastq/$name.$type.R2.fastq.gz > $name/fastq/$name.$type.R2.fastq.gz.md5
	done

	echo '---- Running FastQC before trimming ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	for type in $types;
	do
		fastqc -t $threads \
		$name/fastq/$name.$type.R1.fastq.gz \
		$name/fastq/$name.$type.R2.fastq.gz \
		--outdir=$name/results/QC
	done

	echo '---- Trimming reads ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
	for type in $types;
	do
		trimmomatic_file=$(basename $trimmomatic_dir)
		if [ -z $phred ]; then phred=$(sh $repository_dir/all_DeterminePhred.sh $name $type); fi
		java -Xmx${RAM}G -jar $trimmomatic_dir"/"$trimmomatic_file".jar" PE \
		-threads $threads -$phred \
		$name/fastq/$name.$type.R1.fastq.gz \
		$name/fastq/$name.$type.R2.fastq.gz \
		$temp_dir/$name.$type.R1.passed.fastq.gz \
		$temp_dir/$name.$type.R1.not_passed.fastq.gz \
		$temp_dir/$name.$type.R2.passed.fastq.gz \
		$temp_dir/$name.$type.R2.not_passed.fastq.gz \
		LEADING:25 TRAILING:25 MINLEN:50 \
		SLIDINGWINDOW:10:25 \
		ILLUMINACLIP:$repository_dir/../data/BBDuk-Adapters.fa:2:30:10
	done

	echo '---- Running FastQC after trimming ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	for type in $types;
	do
		fastqc -t $threads \
		$temp_dir/$name.$type.R1.passed.fastq.gz \
		$temp_dir/$name.$type.R2.passed.fastq.gz \
		--outdir=$name/results/QC
	done

	echo '---- Removing fastq files ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
	for type in $types;
	do
		rm $temp_dir/$name.$type.R1.not_passed.fastq.gz
		rm $temp_dir/$name.$type.R2.not_passed.fastq.gz
		rm $name/fastq/$name.$type.R1.fastq.gz
		rm $name/fastq/$name.$type.R2.fastq.gz 
	done

	echo '---- Mapping trimmed reads ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	for type in $types;
	do
		bwa mem -t $threads $genomeindex_dir \
		-Y -K $bwainputbases -v 1 \
		$temp_dir/$name.$type.R1.passed.fastq.gz \
		$temp_dir/$name.$type.R2.passed.fastq.gz \
		| java -Xmx${RAM}G -Dpicard.useLegacyParser=false \
		-jar $picard_dir/picard.jar CleanSam \
		-I /dev/stdin \
		-O $temp_dir/$name.$type.cleaned.bam \
		-VALIDATION_STRINGENCY LENIENT
	done

		for type in $types;
	do
		rm $temp_dir/$name.$type.R1.passed.fastq.gz
		rm $temp_dir/$name.$type.R2.passed.fastq.gz
	done

elif [ $ends = 'SE' ]; then
	echo '---- Calculating md5-sums ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
	for type in $types;
	do
		md5sum $name/fastq/$name.$type.R1.fastq.gz > $name/fastq/$name.$type.R1.fastq.gz.md5
	done

	echo '---- Running FastQC before trimming ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	for type in $types;
	do
		fastqc -t $threads \
		$name/fastq/$name.$type.R1.fastq.gz \
		--outdir=$name/results/QC
	done

	echo '---- Trimming reads ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
	for type in $types;
	do
		trimmomatic_file=$(basename $trimmomatic_dir)
		if [ -z $phred ]; then phred=$(sh $repository_dir/all_DeterminePhred.sh $name $type); fi
		java -Xmx${RAM}G -jar $trimmomatic_dir"/"$trimmomatic_file".jar" SE \
		-threads $threads -$phred \
		$name/fastq/$name.$type.R1.fastq.gz \
		$temp_dir/$name.$type.R1.passed.fastq.gz \
		LEADING:25 TRAILING:25 MINLEN:50 \
		SLIDINGWINDOW:10:25 \
		ILLUMINACLIP:$repository_dir/../data/BBDuk-Adapters.fa:2:30:10
	done

	echo '---- Running FastQC after trimming ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	for type in $types;
	do
		fastqc -t $threads \
		$temp_dir/$name.$type.R1.passed.fastq.gz \
		--outdir=$name/results/QC
	done

	echo '---- Removing fastq files ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
	for type in $types;
	do
		rm $name/fastq/$name.$type.R1.fastq.gz
	done

	echo '---- Mapping trimmed reads ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	for type in $types;
	do
		bwa mem -t $threads $genomeindex_dir \
		-Y -K $bwainputbases -v 1 \
		$temp_dir/$name.$type.R1.passed.fastq.gz \
		| java -Xmx${RAM}G -Dpicard.useLegacyParser=false \
		-jar $picard_dir/picard.jar CleanSam \
		-I /dev/stdin \
		-O $temp_dir/$name.$type.cleaned.bam \
		-VALIDATION_STRINGENCY LENIENT
	done

	for type in $types;
	do
		rm $temp_dir/$name.$type.R1.passed.fastq.gz
	done
fi
	echo '---- Postprocessing I (Sorting, fixing read groups and marking duplicates) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	for type in $types;
	do
		/opt/bin/sambamba sort \
		-t $threads -m ${RAM}GB --tmpdir=$temp_dir \
		-o $temp_dir/$name.$type.cleaned.sorted.bam \
		$temp_dir/$name.$type.cleaned.bam &&

		rm $temp_dir/$name.$type.cleaned.bam &&
		
		java -Xmx${RAM}G -Dpicard.useLegacyParser=false \
		-jar $picard_dir/picard.jar AddOrReplaceReadGroups \
		-I $temp_dir/$name.$type.cleaned.sorted.bam \
		-O $temp_dir/$name.$type.cleaned.sorted.readgroups.bam \
		-ID 1 -LB Lib1 -PL ILLUMINA -PU Run1 -SM $type \
		-MAX_RECORDS_IN_RAM $MAX_RECORDS_IN_RAM &&

		rm $temp_dir/$name.$type.cleaned.sorted.bam &&
		rm $temp_dir/$name.$type.cleaned.sorted.bam.bai &&

		/opt/bin/sambamba markdup \
		--t $threads --tmpdir=$temp_dir \
		--overflow-list-size=$HASH_TABLE_SIZE --hash-table-size=$HASH_TABLE_SIZE \
		$temp_dir/$name.$type.cleaned.sorted.readgroups.bam \
		$temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam &&

		rm $temp_dir/$name.$type.cleaned.sorted.readgroups.bam
	done

	echo '---- Postprocessing II (Base recalibration) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	for type in $types;
	do
		java -Xmx${RAM}G -jar $GATK_dir/gatk.jar BaseRecalibrator \
		-R $genome_file \
		-I $temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam \
		--known-sites $snp_file \
		--use-original-qualities \
		-O $name/results/QC/$name.$type.GATK4.pre.recal.table &&

		java -Xmx${RAM}G -jar $GATK_dir/gatk.jar ApplyBQSR \
		-R $genome_file \
		-I $temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam \
		-O $name/results/bam/$name.$type.bam \
		-bqsr $name/results/QC/$name.$type.GATK4.pre.recal.table &&

		rm $temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam &&
		rm $temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam.bai &&

		java -Xmx${RAM}G -jar $GATK_dir/gatk.jar BaseRecalibrator \
		-R $genome_file \
		-I $name/results/bam/$name.$type.bam \
		--known-sites $snp_file \
		--use-original-qualities \
		-O $name/results/QC/$name.$type.GATK4.post.recal.table &&

		/opt/bin/sambamba index -t $threads $name/results/bam/$name.$type.bam &&

		rm $name/results/bam/$name.$type.bai
	done

if [ $quality_control = "yes" ]; then
	echo '---- Quality control I (Sequencing artifacts, multiple metrics) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	for type in $types;
	do
		java -Xmx${RAM}G -Dpicard.useLegacyParser=false \
		-jar $picard_dir/picard.jar CollectSequencingArtifactMetrics \
		-R $genome_file \
		-I $name/results/bam/$name.$type.bam \
		-O $name/results/QC/$name.$type.bam.artifacts

		java -Xmx${RAM}G -Dpicard.useLegacyParser=false \
		-jar $picard_dir/picard.jar CollectMultipleMetrics \
		-R $genome_file \
		-I $name/results/bam/$name.$type.bam \
		-O $name/results/QC/$name.$type.bam.metrics

		samtools idxstats $name/results/bam/$name.$type.bam \
		> $name/results/QC/$name.$type.bam.idxstats
	done
	
	echo '---- Quality control II (WES- or WGS-specific metrics) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	for type in $types;
	do
		java -Xmx${RAM}G -Dpicard.useLegacyParser=false \
		-jar $picard_dir/picard.jar CollectWgsMetrics \
		-R $genome_file \
		-I $name/results/bam/$name.$type.bam \
		-O $name/results/QC/$name.$type.bam.metrics \
		-SAMPLE_SIZE 100000
	done

	echo '---- Summarizing quality control data ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	multiqc $name/results/QC -n $name -o $name/results/QC/ --pdf --interactive

fi

echo '---- Generate and plot copy number data ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

echo '---- Run CopywriteR ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

Rscript $repository_dir/CNV_RunCopywriter.R $name $species $threads $runmode $genome_dir $centromere_file $varregions_file $types

echo '---- Export raw data and re-normalize using Mode ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

Rscript $repository_dir/CNV_CopywriterGetRawData.R $name $runmode $types
python2 $repository_dir/CNV_CopywriterGetModeCorrectionFactor.py $name
Rscript $repository_dir/CNV_CopywriterGetModeCorrectionFactor.R $name $runmode $types

echo '---- Plot CNV-profiles ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

Rscript $repository_dir/CNV_PlotCopywriter.R $name $species $repository_dir
Rscript $repository_dir/CNV_MapSegmentsToGenes.R $name $species $genecode_file_genes Copywriter 20000 $CGC_file $TruSight_file
sh $repository_dir/CNV_CleanUp.sh $name

echo '---- Run HMMCopy (bin-size 20000) ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

sh $repository_dir/CNV_RunHMMCopy.sh $name $species $config_file $runmode 20000 $types

echo '---- Plot HMMCopy ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

Rscript $repository_dir/CNV_PlotHMMCopy.R $name $species $repository_dir $sequencing_type 20000 \
$mapWig_file $gcWig_file $centromere_file $varregions_file
Rscript $repository_dir/CNV_MapSegmentsToGenes.R $name $species $genecode_file_genes HMMCopy 20000 $CGC_file $TruSight_file

echo '---- Run HMMCopy (bin-size 1000) ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

sh $repository_dir/CNV_RunHMMCopy.sh $name $species $config_file $runmode 1000 $types

echo '---- Finished analysis of sample '$name' ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

exit 0