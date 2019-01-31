#!/bin/bash

##########################################################################################
##
## MouseCancerGenomeAnalysis.sh
##
## Main workflow
##
##########################################################################################

[ $# < 7 ] && { echo "use: sh MouseCancerGenomeAnalysis.sh <name> <fastq_normal_1> <fastq_normal_2> <fastq_tumor_1> <fastq_tumor_2> <sequencing_type> <config_file> [none|CT|GT] [phred33|phred64]"
                 echo "Please edit 'config.sh' before running this script"
                 echo "<name> : Name of the sample"
                 echo "<fastq_normal_1> : path to first normal fastq"
                 echo "<fastq_normal_2> : path to second normal fastq"
                 echo "<fastq_tumor_1> : path to first tumor fastq"
                 echo "<fastq_tumor_2> : path to second tumor fastq"
                 echo "<sequencing_type> : choose from 'WES' or 'WGS'"
                 echo "<config_file> : path to configuration file"
                 echo "<artefact_type> : optional. default none, can be changed to 'GT' (oxidation artefact) or 'CT' (FFPE artefact)"
                 echo "<phred> : optional. default phred33, can be changed to 'phred64'"
                 exit 1
               }

name=$1
fastq_normal_1=$2
fastq_normal_2=$3
fastq_tumor_1=$4
fastq_tumor_2=$5
sequencing_type=$6
config_file=$7
artefact_type=${8:-none}
phred=${9:-phred33}

#reading configuration from $config_file
. $config_file

echo '---- Creating directories ----'
mkdir $temp_dir/
mkdir $name/
mkdir $name/fastq/
mkdir $name/results/
mkdir $name/results/QC
mkdir $name/results/bam
mkdir $name/results/Mutect2
mkdir $name/results/LOH

if [ $sequencing_type = 'WES' ]; then
	mkdir $name/results/Copywriter
elif [ $sequencing_type = 'WGS' ]; then
	mkdir $name/results/Delly
	mkdir $name/results/HMMCopy
	mkdir $name/results/Chromothripsis
fi

echo echo '---- Settings ----' | tee $name/results/QC/$name.report.txt
echo Starting pipeline using these settings: | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt
echo running sample named $name | tee -a $name/results/QC/$name.report.txt
echo using $fastq_normal_1 and $fastq_normal_2 for normal fastqs | tee -a $name/results/QC/$name.report.txt
echo using $fastq_tumor_1 and $fastq_tumor_2 for tumor fastqs | tee -a $name/results/QC/$name.report.txt
echo assuming that experiment is $sequencing_type | tee -a $name/results/QC/$name.report.txt
echo reading configuration file from $config_file | tee -a $name/results/QC/$name.report.txt
echo setting location of repository to $repository_dir | tee -a $name/results/QC/$name.report.txt
echo setting location of genome to $genome_dir | tee -a $name/results/QC/$name.report.txt
echo using mouse reference genome GRCm38.p6 | tee -a $name/results/QC/$name.report.txt
echo setting location for temporary files to $temp_dir| tee -a $name/results/QC/$name.report.txt
echo quality scores are assumed as $phred | tee -a $name/results/QC/$name.report.txt
echo assuming $artefact_type-artefacts for SNV-calling | tee -a $name/results/QC/$name.report.txt
echo starting workflow using $threads Threads and $ram GB of RAM | tee -a $name/results/QC/$name.report.txt

#rerouting STDERR to report file
exec 2>> $name/results/QC/$name.report.txt

echo '---- Copying raw data ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

cp $fastq_normal_1 $name/fastq/$name.Normal.R1.fastq.gz
cp $fastq_normal_2 $name/fastq/$name.Normal.R2.fastq.gz
cp $fastq_tumor_1 $name/fastq/$name.Tumor.R1.fastq.gz
cp $fastq_tumor_2 $name/fastq/$name.Tumor.R2.fastq.gz

echo '---- Running FastQC before trimming ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

fastqc -t $threads \
$name/fastq/$name.Normal.R1.fastq.gz \
$name/fastq/$name.Normal.R2.fastq.gz \
$name/fastq/$name.Tumor.R1.fastq.gz \
$name/fastq/$name.Tumor.R2.fastq.gz \
--outdir=$name/results/QC

echo '---- Trimming reads ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

for type in Normal Tumor;
do
	(
	java -jar $trimmomatic_dir/trimmomatic-0.36.jar PE \
	-threads $threads -$phred \
	$name/fastq/$name.$type.R1.fastq.gz \
	$name/fastq/$name.$type.R2.fastq.gz \
	$temp_dir/$name.$type.R1.passed.fastq \
	$temp_dir/$name.$type.R1.not_passed.fastq \
	$temp_dir/$name.$type.R2.passed.fastq \
	$temp_dir/$name.$type.R2.not_passed.fastq \
	LEADING:25 TRAILING:25 SLIDINGWINDOW:10:25 MINLEN:50
	) &
done

wait

echo '---- Running FastQC after trimming ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

fastqc -t $threads \
$temp_dir/$name.Normal.R1.passed.fastq \
$temp_dir/$name.Normal.R2.passed.fastq \
$temp_dir/$name.Tumor.R1.passed.fastq \
$temp_dir/$name.Tumor.R2.passed.fastq \
--outdir=$name/results/QC

echo '---- Mapping trimmed reads ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

for type in Normal Tumor;
do
	(
	bwa mem -t $threads $genome_dir/bwa_index/GRCm38.p6 \
	$temp_dir/$name.$type.R1.passed.fastq \
	$temp_dir/$name.$type.R2.passed.fastq \
	> $temp_dir/$name.$type.sam
	) &
done

wait

echo '---- Postprocessing I (Cleaning, sorting, fixing read groups and marking duplicates) ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

for type in Normal Tumor;
do
	(
	java -jar $picard_dir/picard.jar CleanSam \
	INPUT=$temp_dir/$name.$type.sam \
	OUTPUT=$temp_dir/$name.$type.cleaned.bam \
	VALIDATION_STRINGENCY=LENIENT

	samtools sort -@ $threads \
	$temp_dir/$name.$type.cleaned.bam \
	-o $temp_dir/$name.$type.cleaned.sorted.bam

	java -jar $picard_dir/picard.jar AddOrReplaceReadGroups \
	I=$temp_dir/$name.$type.cleaned.sorted.bam \
	O=$temp_dir/$name.$type.cleaned.sorted.readgroups.bam \
	ID=1 LB=Lib1-Control PL=ILLUMINA PU=Run1 SM=$type \
	MAX_RECORDS_IN_RAM=50000000

	java -jar $picard_dir/picard.jar MarkDuplicates \
	INPUT=$temp_dir/$name.$type.cleaned.sorted.readgroups.bam \
	OUTPUT=$temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam \
	METRICS_FILE=$name/results/QC/$name.$type.duplicate_metrics.txt \
	REMOVE_DUPLICATES=false ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=50000000
	) &
done

wait

echo '---- Postprocessing II (Base recalibration) ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

for type in Normal Tumor;
do
	(
	java -jar $GATK_dir/gatk.jar BaseRecalibrator \
	-R $genome_dir/GRCm38.p6.fna \
	-I $temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam \
	--known-sites $genome_dir/MGP.v5.snp_and_indels.exclude_wild.vcf.gz \
	--use-original-qualities \
	-O $name/results/QC/$name.$type.GATK4.pre.recal.table

	java -jar $GATK_dir/gatk.jar ApplyBQSR \
	-R $genome_dir/GRCm38.p6.fna \
	-I $temp_dir/$name.$type.cleaned.sorted.readgroups.bam \
	-O $name/results/bam/$name.$type.bam \
	-bqsr $name/results/QC/$name.$type.GATK4.pre.recal.table

	java -jar $GATK_dir/gatk.jar BaseRecalibrator \
	-R $genome_dir/GRCm38.p6.fna \
	-I $name/results/bam/$name.$type.bam \
	--known-sites $genome_dir/MGP.v5.snp_and_indels.exclude_wild.vcf.gz \
	--use-original-qualities \
	-O $name/results/QC/$name.$type.GATK4.post.recal.table

	samtools index -@ $threads $name/results/bam/$name.$type.bam
	) &
done

wait

echo '---- Quality control I (Sequencing artifacts, multiple metrics) ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

for type in Normal Tumor;
do
	(
	java -jar $picard_dir/picard.jar CollectSequencingArtifactMetrics \
	R=$genome_dir/GRCm38.p6.fna \
	I=$name/results/bam/$name.$type.bam \
	O=$name/results/QC/$name.$type.bam.artifacts

	java -jar $picard_dir/picard.jar CollectMultipleMetrics \
	I=$name/results/bam/$name.$type.bam \
	O=$name/results/QC/$name.$type.bam.metrics \
	R=$genome_dir/GRCm38.p6.fna
	) &
done

wait

echo '---- Quality control II (WES- or WGS-specific metrics) ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

if [ $sequencing_type = 'WES' ]; then
	for type in Normal Tumor;
	do
		(
		java -jar $picard_dir/picard.jar CollectHsMetrics \
		SAMPLE_SIZE=100000 \
		I=$name/results/bam/$name.$type.bam \
		O=$name/results/QC/$name.$type.bam.metrics \
		R=$genome_dir/GRCm38.p6.fna \
		BAIT_INTERVALS=$genome_dir/SureSelect_Mouse_All_Exon_V1_mm10.bed.list \
		TARGET_INTERVALS=$genome_dir/SureSelect_Mouse_All_Exon_V1_mm10.bed.list
		) &
	done

	wait

elif [ $sequencing_type = 'WGS' ]; then
	for type in Normal Tumor;
	do
		(
		java -jar $picard_dir/picard.jar CollectWgsMetrics \
		I=$name/results/bam/$name.$type.bam \
		O=$name/results/QC/$name.$type.bam.metrics \
		R=$genome_dir/GRCm38.p6.fna \
		SAMPLE_SIZE=1000000
		) &
	done

	wait
fi

echo '---- Summarizing quality control data ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

multiqc $name/results/QC -n $name -o $name/results/QC/ --pdf

echo '---- Running Mutect2 (matched tumor-normal) ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

java -jar $GATK_dir/gatk.jar Mutect2 \
--native-pair-hmm-threads $threads \
-R $genome_dir/GRCm38.p6.fna \
-I $name/results/bam/$name.Tumor.bam \
-I $name/results/bam/$name.Normal.bam \
-tumor Tumor -normal Normal \
-O $name/results/Mutect2/$name.m2.vcf \
-bamout $name/results/Mutect2/$name.m2.bam

echo '---- Mutect2 Postprocessing I (OrientationFilter, Indel size selection, filtering) ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

java -jar $GATK_dir/gatk.jar FilterMutectCalls \
--variant $name/results/Mutect2/"$name".m2.vcf \
--output $name/results/Mutect2/"$name".m2.filt.vcf

if [ $artefact_type = 'none' ]; then
	cp $name/results/Mutect2/$name.m2.filt.vcf \
	$name/results/Mutect2/$name.m2.filt.AM.vcf
elif [ $artefact_type = 'CT' ]; then
	java -jar $GATK_dir/gatk.jar FilterByOrientationBias \
	-V $name/results/Mutect2/$name.m2.filt.vcf -P \
	$name/results/QC/$name.Tumor.bam.artifacts.pre_adapter_detail_metrics \
	--artifact-modes G/T --output $name/results/Mutect2/$name.m2.filt.AM.vcf
elif [ $artefact_type = 'GT' ]; then
	java -jar $GATK_dir/gatk.jar FilterByOrientationBias \
	-V $name/results/Mutect2/$name.m2.filt.vcf -P \
	$name/results/QC/$name.Tumor.bam.artifacts.pre_adapter_detail_metrics \
	--artifact-modes G/T --output $name/results/Mutect2/$name.m2.filt.AM.vcf
fi

if [ $artefact_type = 'none' ]; then
	cp $name/results/Mutect2/$name.m2.filt.AM.vcf \
	$name/results/Mutect2/$name.m2.filt.AM.filtered.vcf
elif [ $artefact_type = 'CT' ] || [ $artefact_type = 'GT' ]; then
	cat $name/results/Mutect2/$name.m2.filt.AM.vcf \
	| java -jar $snpeff_dir/SnpSift.jar filter \
	"( ( ( na FILTER ) & (exists GEN[Tumor].OBP) & \
	(GEN[Tumor].OBP <= 0.05) ) | ( ( na FILTER ) ) )" \
	> $name/results/Mutect2/$name.m2.filt.AM.filtered.vcf
fi

java -jar $GATK_dir/gatk.jar SelectVariants --max-indel-size 10 \
-V $name/results/Mutect2/$name.m2.filt.AM.filtered.vcf \
-output $name/results/Mutect2/$name.m2.filt.AM.filtered.selected.vcf

cat $name/results/Mutect2/$name.m2.filt.AM.filtered.selected.vcf \
| java -jar $snpeff_dir/SnpSift.jar filter \
"( ( na FILTER ) & (GEN[Tumor].AF >= 0.1) & \
( ( GEN[Tumor].AD[0] + GEN[Tumor].AD[1]) >= 10 ) & \
( ( GEN[Normal].AD[0] + GEN[Normal].AD[1]) >= 10 ) & \
(GEN[Tumor].AD[1] > 2) & (GEN[Normal].AD[1] = 0) )" \
> $name/results/Mutect2/$name.m2.postprocessed.vcf

echo '---- Mutect2 Postprocessing II (Filtering out known SNV using Sanger database) ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

bgzip $name/results/Mutect2/$name.m2.postprocessed.vcf

tabix -p vcf $name/results/Mutect2/$name.m2.postprocessed.vcf.gz

bcftools isec -C -c none -O z -w 1 \
-o $name/results/Mutect2/$name.m2.postprocessed.snp_removed.vcf.gz \
$name/results/Mutect2/$name.m2.postprocessed.vcf.gz \
$genome_dir/MGP.v5.snp_and_indels.exclude_wild.vcf.gz

mv $name/results/Mutect2/$name.m2.postprocessed.snp_removed.vcf.gz \
$name/results/Mutect2/$name.Mutect2.vcf.gz

gunzip $name/results/Mutect2/$name.Mutect2.vcf.gz

echo '---- Mutect2 Postprocessing III (Annotate calls) ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

java -Xmx"$ram"g -jar $snpeff_dir/snpEff.jar GRCm38.86 -canon -v \
-csvStats $name/results/Mutect2/$name.Mutect2.annotated.vcf.stats \
$name/results/Mutect2/$name.Mutect2.vcf \
> $name/results/Mutect2/$name.Mutect2.annotated.vcf

cat $name/results/Mutect2/$name.Mutect2.annotated.vcf \
| $snpeff_dir/scripts/vcfEffOnePerLine.pl \
> $name/results/Mutect2/$name.Mutect2.annotated.one.vcf

java -jar $snpeff_dir/SnpSift.jar extractFields \
$name/results/Mutect2/$name.Mutect2.annotated.one.vcf \
CHROM POS REF ALT "GEN[Tumor].AF" "GEN[Tumor].AD[0]" "GEN[Tumor].AD[1]" \
"GEN[Normal].AD[0]" "GEN[Normal].AD[1]" ANN[*].GENE  ANN[*].EFFECT \
ANN[*].IMPACT ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
> $name/results/Mutect2/$name.Mutect2.txt

sh $repository_dir/SNV_CleanUp.sh $name MS

echo '---- Running Mutect2 (single-sample tumor and normal) ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

for type in Normal Tumor;
do
	(
	java -jar $GATK_dir/gatk.jar Mutect2 \
	--native-pair-hmm-threads $threads \
	-R $genome_dir/GRCm38.p6.fna \
	-I $name/results/bam/$name.$type.bam \
	-tumor $type \
	-O $name/results/Mutect2/$name."$type".m2.vcf
	-bamout $name/results/Mutect2/$name."$type".m2.bam
	) &
done

wait

echo '---- Mutect2 Postprocessing for LOH ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

for type in Normal Tumor;
do
	(
	java -jar $GATK_dir/gatk.jar FilterMutectCalls \
	--variant $name/results/Mutect2/"$name".$type.m2.vcf \
	--output $name/results/Mutect2/"$name".$type.m2.filt.vcf

	java -jar $snpeff_dir/SnpSift.jar extractFields \
	$name/results/Mutect2/"$name".$type.m2.filt.vcf \
	CHROM POS REF ALT "GEN["$type"].AF" "GEN["$type"].AD[0]" \
	"GEN["$type"].AD[1]" "GEN["$type"].MMQ" "GEN["$type"].MBQ" \
	"GEN["$type"].SA_POST_PROB[0]" "GEN["$type"].SA_POST_PROB[1]" \
	"GEN["$type"].SA_POST_PROB[2]" \
	> $name/results/Mutect2/"$name"."$type".Mutect2.Positions.txt
	) &
done

wait

sh $repository_dir/SNV_CleanUp.sh $name SS

echo '---- Generate LOH data ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

Rscript $repository_dir/LOH_GenerateVariantTable.R \
$name $genome_dir/GRCm38.p6.fna $repository_dir

Rscript $repository_dir/LOH_MakePlots.R \
$name $repository_dir

echo '---- Generate and plot copy number data ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

if [ $sequencing_type = 'WES' ]; then
	echo '---- Run CopywriteR ----' | tee -a $name/results/QC/$name.report.txt
	date | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/CNV_RunCopywriter.R $name $threads

	echo '---- Export raw data and re-normalize using Mode ----' | tee -a $name/results/QC/$name.report.txt
	date | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/CNV_CopywriterGetRawData.R $name
	python $repository_dir/CNV_CopywriterGetModeCorrectionFactor.py $name
	Rscript $repository_dir/CNV_CopywriterGetModeCorrectionFactor.R $name

	echo '---- Plot CNV-profiles ----' | tee -a $name/results/QC/$name.report.txt
	date | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/CNV_PlotCopywriter.R $name $repository_dir
	sh $repository_dir/CNV_CleanUp $name

elif [ $sequencing_type = 'WGS' ]; then
	echo '---- Run HMMCopy ----' | tee -a $name/results/QC/$name.report.txt
	date | tee -a $name/results/QC/$name.report.txt

	sh $repository_dir/CNV_RunHMMCopy.sh $name 20000
	Rscript $repository_dir/CNV_PlotHMMCopy.R $name $repository_dir 20000
fi

echo '---- Optional for WGS: Run Delly ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

if [ $sequencing_type = 'WGS' ]; then
	$delly_dir/delly_v0.7.8_parallel_linux_x86_64bit call \
	-o $name/results/Delly/$name.pre.bcf \
	-g $genome_dir/GRCm38.p6.fna \
	$name/results/bam/$name.Tumor.bam \
	$name/results/bam/$name.Normal.bam

	$delly_dir/delly_v0.7.8_parallel_linux_x86_64bit filter \
	-f somatic -o $name/results/Delly/$name.bcf \
	-s $genome_dir/Samples.tsv $name/results/Delly/$name.pre.bcf

	bcftools view $name/results/Delly/$name.pre.bcf \
	> $name/results/Delly/$name.pre.vcf
fi

echo '---- Optional for WGS: Infer chromothripsis ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

format="tif"

echo '---- Preparing input files and calculating coverage ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

coverage=$(sh $repository_dir/Chromothripsis_GetCoverage.sh $name)

sh $repository_dir/Chromothripsis_FormatTable.sh $name

Rscript $repository_dir/Chromothripsis_AnnotateRatios.R \
-i $name/results/Delly/$name.breakpoints.tab \
 > $name/results/Delly/$name.breakpoints_annotated.tab

Rscript $repository_dir/Chromothripsis_FilterDelly.R \
-n $name -c $coverage \
-i $name/results/Delly/$name.breakpoints_annotated.tab

for chr in $( seq 19 ); do
	if [ $(Rscript $repository_dir/Chromothripsis_RearrangementCounter.R -i $name/results/Delly/$name.breakpoints.filtered.tab -c $chr) -ge 4 ]; then
		echo 'Analysing Chromosome '$chr
		mkdir $name'/results/Chromothripsis/Chr'$chr
		echo '---- Hallmark: Clustering of breakpoints for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		date | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_DetectBreakpointClustering.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-c $chr -n $name -f $format

		echo '---- Hallmark: Regularity of oscillating copy number states for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		date | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_SimulateCopyNumberStates.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-o mouse -c $chr -n $name -s 1000 -a 1000 -f $format -v 1

		echo '---- Hallmark: Interspersed loss and retention of heterozygosity for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		date | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_PlotLOHPattern.R \
		-s $name/results/HMMCopy/$name.HMMCopy.$resolution.segments.txt \
		-d $name/results/HMMCopy/$name.HMMCopy.$resolution.log2RR.txt \
		-v $name/results/LOH/$name.VariantsForLOH.txt \
		-o mouse -c $chr -n $name -f $format

		echo '---- Hallmark: Randomness of DNA fragment joins and segment order for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		date | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_DetectRandomJoins.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-c $chr -n $name -f $format

		echo '---- Hallmark: Ability to walk the derivative chromosome for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		date | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_WalkDerivativeChromosome.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-c $chr -n $name -f $format

		echo '---- Visualisation: Copy number profile combined with complex rearrangements for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		date | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_PlotRearrangementGraph.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-d $name/results/HMMCopy/$name.HMMCopy.$resolution.log2RR.txt \
		-c $chr -n $name -f $format
	else
		echo 'There are too few rearrangements in chromosome '$chr'.'
	fi
done

echo '---- Removing temporary files ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

rm -r $temp_dir
rm $name/results/bam/$name.Normal.bai $name/results/bam/$name.Tumor.bai

echo '---- Run msisensor----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

mkdir $name/results/msisensor
msisensor msi -n $name/results/bam/$name.Normal.bam \
-t $name/results/bam/$name.Tumor.bam \
-o $name/results/msisensor/"$name".msisensor \
-d $genome_dir/GRCm38.p6.microsatellites

echo '---- "Remove raw fastq files (Y/N, followed by ENTER)?" ----' | tee -a $name/results/QC/$name.report.txt
read -p "Remove raw fastq files (Y/N)?" remove_raw_files
case $remove_raw_files in
Y|y) echo '---- Removing raw sequencing files ----' | tee -a $name/results/QC/$name.report.txt;
	 date | tee -a $name/results/QC/$name.report.txt;
     rm -r $name/fastq;;
N|n) exit;;
esac

echo '---- Finished analysis of sample $name ----' | tee -a $name/results/QC/$name.report.txt
date | tee -a $name/results/QC/$name.report.txt

exit 0