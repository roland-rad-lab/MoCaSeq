#!/bin/bash

##########################################################################################
##
## MoCaSeq.sh
##
## Main workflow
##
##########################################################################################

usage()
{
	echo "  Usage: $0 "
	echo "	-n, --name               Name of the sample."
	echo "	-s, --species            Set to 'Mouse' or 'Human'. Defaults to Mouse."
	echo "	-nf, --fastq_normal_fw   Path to first normal Fastq. Do NOT use if running single-sample tumor only."
	echo "	-nr, --fastq_normal_rev  Path to second normal Fastq. Do NOT use if running single-sample tumor only."
	echo "	-tf, --fastq_tumor_fw    Path to first tumor fastq. Do NOT use if running single-sample normal only."
	echo "	-tr, --fastq_tumor_rev   Path to second tumor fastq. Do NOT use if running single-sample normal only."
	echo "	-nb, --bam_normal        Path to normal BAM. Do NOT use in combination with -nf or -nr. When used, -rm MUST be specified."
	echo "	-tb, --bam_tumor         Path to tumor BAM. Do NOT use in combination with -tf or -tr. When used, -rm MUST be specified."
	echo "	-rm, --repeat_mapping    If -nb or -tb are specified, determines whether mapping is re-done ('yes') or whether the complete mapping procedure is skipped ('no')."
	echo "	-st, --sequencing_type   Set to 'WES' or 'WGS'. Defaults to WES."
	echo "	-c, --config             Path to configuration file. Optional."
	echo "	-qc, --quality_control   Determines whether QC is done ('yes') or skipped ('no'). Optional."
	echo "	-t, --threads            Number of CPU threads. Optional. Defaults to 8."
	echo "	-r, --RAM                Amount of Gb RAM. Optional. Defaults to 32."
	echo "	-temp, --temp_dir        Path to temporary directory. Optional. Defaults to current working directory."
	echo "	-art, --artefact         Set to 'no' for no filter or 'yes' to automatically filter read-orientation bias artifacts using GATK (ob-priors). This includes OxoG oxidation artefacts and FFPE artefacts. Optional. Defaults to yes."
	echo "	-filt, --filtering       Set to 'soft' (AF >= 0.05, , Variant in Tumor >= 2, Variant in Normal <= 1, Coverage >= 5, dbSNP common for human), 'hard' (AF >= 0.1, Variant in Tumor >= 3, Variant in Normal = 0, Coverage >= 10, dbSNP all for human) or 'none' (no filters). Optional. Defaults to 'soft'."
	echo "	-p, --phred              If not set, script will try to automatically extract phred-score. Otherwise, set manually to 'phred33' or 'phred64'. 'phred64' only relevant for Illumina data originating before 2011. Optional."
	echo "	-mu, --Mutect2           Set to 'yes' or 'no'. Needed for LOH analysis and Titan. Greatly increases runtime for WGS. Optional. Defaults to 'yes'."
	echo "	-ck, --CNVKit            Set to 'yes' or 'no'. Needed for LOH analysis. Optional. Defaults to 'yes'."
	echo "	-de, --Delly             Set to 'yes' or 'no'. Needed for chromothripsis inference. Do not use for WES. Optional. Defaults to 'no'. Only use in matched-sample mode."
	echo "	-ti, --Titan             Set to 'yes' or 'no'. Runs TITAN to model subclonal copy number alterations, predict LOH and estimate tumor purity. Greatly increases runtime for WGS. If set to 'yes', forces Mutect2 to 'yes'. Optional. Defaults to 'yes' for WES and 'no' for WGS. Only use in matched-sample mode."
	echo "	-abs, --Absolute         Set to 'yes' or 'no'. Runs ABSOLUTE to estimate purity/ploidy and compute copy-numbers. Optional. Can also include information from somatic mutation data, for this set Mutect2 to 'yes'. Defaults to 'yes' for WES and WGS."
	echo "	-fac, --Facets           Set to 'yes' or 'no'. Runs the allele-specific copy-number caller FACETS with sample purity estimations. Optional. Defaults to 'yes' for WES and 'no' for WGS. Only use in matched-sample mode."
	echo "	-bt, --BubbleTree             Set to 'yes' or 'no'. Runs the analysis of tumoral aneuploidy and clonality. Optional. If set to 'yes', forces Mutect2 to 'yes'. Optional. Defaults to 'yes' for WES and WGS. Only use in matched-sample mode."
	echo "	-gatk, --GATKVersion     Set to '4.1.7.0' or '4.2.0.0', determining which GATK version is used. Optional. Defaults to 4.2.0.0 (high recommended to identify all mutations)"
	echo "	--test                   If set to 'yes': Will download reference files (if needed) and start a test run. All other parameters will be ignored"
	echo "	--memstats               If integer > 0 specified, will write timestamped memory usage and cumulative CPU time usage of the docker container to ./results/memstats.txt every <integer> seconds. Defaults to '0'."
	echo "  --para                   Run Mutect2 in parallel"
	echo "	--help                   Show this help."
  exit 1
}

# quick and dirty charliecloud solutions
cd /var/pipeline
PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/bin

# default parameters
fastq_normal_1=
fastq_normal_2=
fastq_tumor_1=
fastq_tumor_2=
bam_normal=
bam_tumor=
repeat_mapping=yes
sequencing_type=WES
quality_control=yes
threads=8
RAM=32
temp_dir=/var/pipeline/temp
artefact_type=yes
filtering=soft
phred=
Mutect2=yes
CNVKit=yes
Delly=no
runmode=MS
GATK=4.2.0.0
test=no
memstats=0
config_file=
species=Mouse
para=no
# Titan and Facets will be set below
BubbleTree=yes
Absolute=yes

# parse parameters
if [ "$1" = "" ]; then usage; fi
while [ "$1" != "" ]; do case $1 in
	-n|--name) shift;name="$1";;
	-s|--species) shift;species="$1";;
	-nf|--fastq_normal_fw) shift;fastq_normal_1="$1";;
	-nr|--fastq_normal_rev) shift;fastq_normal_2="$1";;
	-tf|--fastq_tumor_fw) shift;fastq_tumor_1="$1";;
	-tr|--fastq_tumor_rev) shift;fastq_tumor_2="$1";;
	-nb|--bam_normal) shift;bam_normal="$1";;
	-tb|--bam_tumor) shift;bam_tumor="$1";;
	-rm|--repeat_mapping) shift;repeat_mapping="$1";;
	-rq|--quality_control) shift;quality_control="$1";;
	-st|--sequencing_type) shift;sequencing_type="$1";;
	-c|--config) shift;config_file="$1";;
	-t|--threads) shift;threads="$1";;
	-r|--RAM) shift;RAM="$1";;
	-temp|--temp_dir) shift;temp_dir="$1";;
	-art|--artefact) shift;artefact_type="$1";;
	-filt|--filtering) shift;filtering="$1";;
    -p|--phred) shift;phred="$1";;
    -mu|--Mutect2) shift;Mutect2="$1";;
    -ck|--CNVKit) shift;CNVKit="$1";;
    -de|--Delly) shift;Delly="$1";;
    -ti|--Titan) shift;Titan="$1";;
		-abs|--Absolute) shift;Absolute="$1";;
		-fac|--Facets) shift;Facets="$1";;
		-bt|--BubbleTree) shift;BubbleTree="$1";;
    -gatk|--GATKVersion) shift;GATK="$1";;
    --memstats) shift;memstats="$1";;
    --test) shift;test="$1";;
		--para) shift;para="$1";;
    --help) usage;shift;;
	*) usage;shift;;
esac; shift; done

if [ -z $config_file ]; then
	config_file=/opt/MoCaSeq/config.sh
fi

# init variables of test run
test_dir=${config_file%/*}/test

# function to check if a file exists and is not empty
function CheckFile {
  file=$1

  if [ ! -s "$file" ]
  then
   echo "ERROR! File not created successfully: $file"
   exit 1
   fi
}


#test=yes
if [ $test = 'yes' ]; then
	name=MoCaSeq_Test
	species=Mouse
	fastq_normal_1=$test_dir/Mouse.Normal.R1.fastq.gz
	fastq_normal_2=$test_dir/Mouse.Normal.R2.fastq.gz
	fastq_tumor_1=$test_dir/Mouse.Tumor.R1.fastq.gz
	fastq_tumor_2=$test_dir/Mouse.Tumor.R2.fastq.gz
	sequencing_type=WES
	bam_normal=
	bam_tumor=
	repeat_mapping=yes
	quality_control=yes
	threads=4
	RAM=8
	temp_dir=/var/pipeline/temp
	artefact_type=no
	filtering=hard
	phred=
	Mutect2=yes
	CNVKit=no
	Titan=no
	Absolute=no
	Facets=no
	BubbleTree=no
	Delly=no
	runmode=MS
fi

# set some species specific arguments
if [ $species = 'Mouse' ]; then
	echo 'Species set to Mouse'
	chromosomes=19
	echo "Species $species_lowercase"
elif [ $species = 'Human' ]; then
	echo 'Species set to Human'
	chromosomes=22
else echo "Invalid species input (${species}). Choose Mouse or Human"; exit 1
fi

# TODO: add parameter check (in general, for all values)
if [ $filtering = 'no' ]; then
	filtering="none"
fi

species_lowercase=${species,,} # set the species to lowercase to match some scripts input format (e.g. Chromothripsis)

if [ -z $fastq_normal_1 ] && [ -z $fastq_normal_2 ] && [ ! -z $fastq_tumor_1 ] && [ ! -z $fastq_tumor_2 ] && [ -z $bam_normal ] && [ -z $bam_tumor ]; then
	runmode="SS"
	types="Tumor"
elif [ ! -z $fastq_normal_1 ] && [ ! -z $fastq_normal_2 ] && [ -z $fastq_tumor_1 ] && [ -z $fastq_tumor_2 ] && [ -z $bam_normal ] && [ -z $bam_tumor ]; then
	runmode="SS"
	types="Normal"
elif [ ! -z $fastq_normal_1 ] && [ ! -z $fastq_normal_2 ] && [ ! -z $fastq_tumor_1 ] && [ ! -z $fastq_tumor_2 ] && [ -z $bam_normal ] && [ -z $bam_tumor ]; then
	runmode="MS"
	types="Tumor Normal"
elif [ -z $fastq_normal_1 ] && [ -z $fastq_normal_2 ] && [ -z $fastq_tumor_1 ] && [ -z $fastq_tumor_2 ] && [ -z $bam_normal ] && [ ! -z $bam_tumor ]; then
	runmode="SS"
	types="Tumor"
elif [ -z $fastq_normal_1 ] && [ -z $fastq_normal_2 ] && [ -z $fastq_tumor_1 ] && [ -z $fastq_tumor_2 ] && [ ! -z $bam_normal ] && [ -z $bam_tumor ]; then
	runmode="SS"
	types="Normal"
elif [ -z $fastq_normal_1 ] && [ -z $fastq_normal_2 ] && [ -z $fastq_tumor_1 ] && [ -z $fastq_tumor_2 ] && [ ! -z $bam_normal ] && [ ! -z $bam_tumor ]; then
	runmode="MS"
	types="Tumor Normal"
else echo 'Invalid combination of input files. Either use -tf/-tr/-nf/-nr OR -tb/-nb'; #exit 1
fi

# # CHECK IF ALL NEEDED FILE EXIST
# if [ "$fastq_normal_1" != "" ] && [ ! -f "$fastq_normal_1" ]; then
# echo "ERROR, File not found: $fastq_normal_1"
# exit 1
# fi
#
# if [ "$fastq_normal_2" != "" ] && [ ! -f "$fastq_normal_2" ]; then
# echo "ERROR, File not found: $fastq_normal_2"
# exit 1
# fi
#
# if [ "$fastq_tumor_1" != "" ] && [ ! -f "$fastq_tumor_1" ]; then
# echo "ERROR, File not found: $fastq_tumor_1"
# exit 1
# fi
#
# if [ "$fastq_tumor_2" != "" ] && [ ! -f "$fastq_tumor_2" ]; then
# echo "ERROR, File not found: $fastq_tumor_2"
# exit 1
# fi
#
# if [ "$bam_normal" != "" ] && [ ! -f "$bam_normal" ]; then
# echo "ERROR, File not found: $bam_normal"
# exit 1
# fi
#
# if [ "$bam_tumor" != "" ] && [ ! -f "$bam_tumor" ]; then
# echo "ERROR, File not found: $bam_tumor"
# exit 1
# fi

if [ "$GATK" == "4.1.0.0" ]; then
echo "ERROR, GATK 4.1.0.0 is not supported anymore, please choose 4.2.0.0"
exit 1
fi


# SET PARAMETERS FOR PURITY ANALYSIS

# TITAN ist default 'yes' for WES and default 'no' for WGS
if [ $sequencing_type = 'WES' ] && [ -z $Titan ]; then
	Titan=yes
elif [ $sequencing_type = 'WGS' ] && [ -z $Titan ]; then
	Titan=no
fi

# If TITAN set to 'yes', forces Mutect2 to 'yes'. set to 'no' if runmode SS (needs tumor+normal WIG files)
if [ $Titan = 'yes' ] && [ $runmode = 'MS' ]; then
	Mutect2=yes
	echo 'PARAMETER CHANGE: TITAN needs LOH data from Mutect2. Mutect2 was set to "yes".'
	Titan=yes
else
	if [ $Titan = 'yes' ]; then
	echo 'PARAMETER CHANGE: TITAN needs matched tumor/normal and can not be used on single sample runs. TITAN was set to "no".'
	fi
	Titan=no
fi


# Facets ist default 'yes' for WES and default 'no' for WGS
if [ $sequencing_type = 'WES' ] && [ -z $Facets ]; then
	Facets=yes
elif [ $sequencing_type = 'WGS' ] && [ -z $Facets ]; then
	Facets=no
fi

# Facets needs tumor and normal BAM files, so it can only be used in MS mode
if [ $Facets = 'yes' ] && [ $runmode = 'SS' ]; then
	Facets=no
	echo 'PARAMETER CHANGE: FACETS needs matched tumor/normal and can not be used on single sample runs. Facets was set to "no".'
fi


# BubbleTree needs LOH and paired samples (tumor+normal WIG files)
# If BubbleTree set to 'yes', forces Mutect2 to 'yes'. set to 'no' if runmode SS
if [ $BubbleTree = 'yes' ] && [ $runmode = 'MS' ]; then

	if [ $Mutect2 = 'no' ]; then
	Mutect2=yes
	echo 'PARAMETER CHANGE: BubbleTree needs LOH data from Mutect2. Mutect2 was set to "yes".'
	fi

	BubbleTree=yes
else

	if [ $BubbleTree = 'yes' ]; then
	echo 'PARAMETER CHANGE: BubbleTree needs matched tumor/normal and can not be used on single sample runs. BubbleTree was set to "no".'
	fi

	BubbleTree=no
fi

# check for valid entries for the bool variables
boolVars="repeat_mapping quality_control artefact_type Mutect2 CNVKit Delly BubbleTree Facets Titan Absolute test"
for varName in $boolVars
do
	eval varValue=\$$varName
	if [ $varValue != 'yes' ] && [ $varValue != 'no' ]; then
		echo "Invalid value for ${varName} (only yes or no are allowed): ${varValue}"
	fi
done

if [ $filtering != 'soft' ] && [ $filtering != 'hard' ]; then
	echo "Invalid value for filtering (only soft or hard are allowed): ${filtering}"
fi


#reading configuration from $config_file
source $config_file
repository_dir=${config_file%/*}/repository

echo "---- Starting ${species} Cancer Genome Analysis ----"
echo -e "$(date) \t timestamp: $(date +%s)"

echo '---- Creating directories ----'
echo -e "$(date) \t timestamp: $(date +%s)"
mkdir -p $name/results/QC
mkdir -p $name/results/bam
mkdir -p $name/fastq


# log memory and cpu usage
logstats(){
	echo -e "date \t timestamp \t memory_usage_bytes \t cumulative_cpu_nanoseconds \t cores" > $name/results/memstats.txt
	while sleep $memstats; do ($repository_dir/Meta_logstats.sh >> $name/results/memstats.txt &) ; done
}

if [ $memstats -gt 0 ]; then
	logstats &
fi


echo '---- Checking for available reference files ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

# Niklas: removed the "-z ...txt" since it should be "! -f" to work and therefore was not used anyways
if grep --quiet DONE $genome_dir/GetReferenceData.txt; then
	echo '---- Reference files not found - Files will be downloaded ----' | tee -a $name/results/QC/${name}.report.txt
	# echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/${name}.report.txt
  # rm -rf $genome_dir
	# if [ $species = 'Mouse' ]; then
	# sh $repository_dir/Preparation_GetReferenceDataMouse.sh $config_file $temp_dir
	# elif [ $species = 'Human' ]; then
	# sh $repository_dir/Preparation_GetReferenceDataHuman.sh $config_file $temp_dir
	# fi
else
	echo '---- Reference files found! ----' | tee -a ${name}/results/QC/${name}.report.txt
fi

# check if all needed files are given in the reference folder
# || exit 1 will exit if the secondary script fails and itself calls "exit 1"
$repository_dir/CheckReferenceFiles.sh $FileList  || exit 1

# CREATE SUBDIRS
if [ ! -d $temp_dir ]; then
  mkdir -p $temp_dir/
fi

if [ $RAM -ge 16 ]; then
	bwainputbases=100000000
else bwainputbases=10000000
fi

MAX_RECORDS_IN_RAM=$(expr $RAM \* 250000)
HASH_TABLE_SIZE=$((RAM*1000000000/500))

echo "---- Starting ${species} Cancer Genome Analysis ----" | tee -a $name/results/QC/$name.report.txt
echo Starting pipeline using these settings: | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
echo Running sample named $name | tee -a $name/results/QC/$name.report.txt
echo Running in $runmode-mode | tee -a $name/results/QC/$name.report.txt

if [ $runmode = "MS" ] && [ ! -z $fastq_normal_1 ] && [ ! -z $fastq_normal_2 ] && [ ! -z $fastq_tumor_1 ] && [ ! -z $fastq_tumor_2 ] && [ -z $bam_normal ] && [ -z $bam_tumor ]; then
	echo Using $fastq_normal_1 and $fastq_normal_2 for normal fastqs | tee -a $name/results/QC/$name.report.txt
	echo Using $fastq_tumor_1 and $fastq_tumor_2 for tumor fastqs | tee -a $name/results/QC/$name.report.txt
elif [ $runmode = "MS" ] && [ -z $fastq_normal_1 ] && [ -z $fastq_normal_2 ] && [ -z $fastq_tumor_1 ] && [ -z $fastq_tumor_2 ] && [ ! -z $bam_normal ] && [ ! -z $bam_tumor ]; then
	echo Using $bam_normal for normal bam | tee -a $name/results/QC/$name.report.txt
	echo Using $bam_tumor for tumor bam | tee -a $name/results/QC/$name.report.txt
elif [ $runmode = "SS" ] && [ $repeat_mapping = "no" ] && [ -z $fastq_normal_1 ] && [ -z $fastq_normal_2 ]  && [ -z $bam_normal ] && [ -z $bam_tumor ]; then
	echo Assigning $fastq_tumor_1 and $fastq_tumor_2 as $types | tee -a $name/results/QC/$name.report.txt
elif [ $runmode = "SS" ] && [ $repeat_mapping = "no" ] && [ -z $fastq_tumor_1 ] && [ -z $fastq_tumor_2 ]  && [ -z $bam_normal ] && [ -z $bam_tumor ]; then
	echo Assigning $fastq_normal_1 and $fastq_normal_2 as $types | tee -a $name/results/QC/$name.report.txt
elif [ $runmode = "SS" ] && [ $repeat_mapping = "yes" ] && [ -z $bam_normal ] ; then
	echo Assigning $bam_tumor as $types | tee -a $name/results/QC/$name.report.txt
elif [ $runmode = "SS" ] && [ $repeat_mapping = "yes" ] && [ -z $bam_tumor ] ; then
	echo Assigning $bam_normal as $types | tee -a $name/results/QC/$name.report.txt
fi
if [ $repeat_mapping = "no" ]; then
	echo Input BAMs will NOT be re-mapped | tee -a $name/results/QC/$name.report.txt
fi
echo Assuming that reads are from $species | tee -a $name/results/QC/$name.report.txt
echo Assuming that experiment is $sequencing_type | tee -a $name/results/QC/$name.report.txt
echo Reading configuration file from $config_file | tee -a $name/results/QC/$name.report.txt
echo Setting location of repository to $repository_dir | tee -a $name/results/QC/$name.report.txt
echo Setting location of genome to $genome_dir | tee -a $name/results/QC/$name.report.txt
echo Setting location for temporary files to $temp_dir| tee -a $name/results/QC/$name.report.txt
echo Filtering orientation bias artefacts for SNV-calling: $artefact_type | tee -a $name/results/QC/$name.report.txt
echo $filtering is setting for filtering of SNV calls | tee -a $name/results/QC/$name.report.txt

echo Quality scores are assumed as $phred | tee -a $name/results/QC/$name.report.txt



#rerouting STDERR to report file
exec 2>> $name/results/QC/$name.report.txt
# exec 2>&1 | tee $name/results/QC/$name.report.txt # print to log file and console





echo '---- Creating directories ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

echo '---- Copying raw data ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

checkR1fq=$name/fastq/${name}.${types}.R1.fastq.gz
checkR2fq=$name/fastq/${name}.${types}.R2.fastq.gz
if [ -f "${checkR1fq}" ] && [ -f "${checkR2fq}" ]; then
	echo -e "skipping sam2fasq step, found pre-existing files:\n ${checkR1fq}\n${checkR2fq}"
else
	# 1. this will either copy fastqs to location OR do sam2fastq
	if [ -z $bam_normal ] && [ -z $bam_tumor ]; then

		echo '---- Copying FASTQs ----' | tee -a $name/results/QC/$name.report.txt

		if [ $runmode = 'MS' ]; then
		cp $fastq_normal_1 $name/fastq/$name.Normal.R1.fastq.gz
		cp $fastq_normal_2 $name/fastq/$name.Normal.R2.fastq.gz
		cp $fastq_tumor_1 $name/fastq/$name.Tumor.R1.fastq.gz
		cp $fastq_tumor_2 $name/fastq/$name.Tumor.R2.fastq.gz

		elif [ $runmode = 'SS' ] && [ $types = 'Tumor' ]; then
		cp $fastq_tumor_1 $name/fastq/$name.$types.R1.fastq.gz
		cp $fastq_tumor_2 $name/fastq/$name.$types.R2.fastq.gz

		elif [ $runmode = 'SS' ] && [ $types = 'Normal' ]; then
		cp $fastq_normal_1 $name/fastq/$name.$types.R1.fastq.gz
		cp $fastq_normal_2 $name/fastq/$name.$types.R2.fastq.gz
		fi

	elif [ $repeat_mapping = "yes" ]; then
		echo '     ---- Converting Bam to Fastq ----' | tee -a $name/results/QC/$name.report.txt

		# remove unmapped reads (to avoid "Mapq Should Be 0 For Unmapped Read")
		#samtools view -bF 4 /var/fastqs/EGAF00001721862/PCSI_0612_Ag_M_526.bam > filtered.bam
		#samtools view -bF 4 /var/fastqs/EGAF00001721862/PCSI_0612_Ag_M_526.bam > filtered.bam

		# BAM files have to be sorted! We will assume they are (because of runtime) but this would be the code to fix it:
		# samtools sort $bam_tumor > $temp_dir/${name}.Tumor.raw.sorted.bam
		# samtools index $temp_dir/${name}.Tumor.raw.sorted.bam
		# bam_tumor=${name}.Tumor.sorted.bam
		# samtools sort $bam_normal > $temp_dir/${name}.Normal.raw.sorted.bam
		# samtools index $temp_dir/${name}.Normal.raw.sorted.bam
		# bam_normal=${name}.Normal.sorted.bam

		if [ $runmode = 'MS' ]; then
		java -Xmx${RAM}G -Dpicard.useLegacyParser=false -jar $picard_dir/picard.jar SamToFastq \
		-INPUT $bam_tumor \
		-FASTQ $name/fastq/$name.Tumor.R1.fastq.gz \
		-SECOND_END_FASTQ $name/fastq/$name.Tumor.R2.fastq.gz \
		-INCLUDE_NON_PF_READS true \
		-VALIDATION_STRINGENCY LENIENT

		java -Xmx${RAM}G -Dpicard.useLegacyParser=false -jar $picard_dir/picard.jar SamToFastq \
		-INPUT $bam_normal \
		-FASTQ $name/fastq/$name.Normal.R1.fastq.gz \
		-SECOND_END_FASTQ $name/fastq/$name.Normal.R2.fastq.gz \
		-INCLUDE_NON_PF_READS true \
		-VALIDATION_STRINGENCY LENIENT

		elif [ $runmode = 'SS' ] && [ $types = 'Tumor' ]; then
		java -Xmx${RAM}G -Dpicard.useLegacyParser=false -jar $picard_dir/picard.jar SamToFastq \
		-INPUT $bam_tumor \
		-FASTQ $name/fastq/$name.Tumor.R1.fastq.gz \
		-SECOND_END_FASTQ $name/fastq/$name.Tumor.R2.fastq.gz \
		-INCLUDE_NON_PF_READS true \
		-VALIDATION_STRINGENCY LENIENT

		elif [ $runmode = 'SS' ] && [ $types = 'Normal' ]; then
		java -Xmx${RAM}G -Dpicard.useLegacyParser=false -jar $picard_dir/picard.jar SamToFastq \
		-INPUT $bam_normal \
		-FASTQ $name/fastq/$name.Normal.R1.fastq.gz \
		-SECOND_END_FASTQ $name/fastq/$name.Normal.R2.fastq.gz \
		-INCLUDE_NON_PF_READS true \
		-VALIDATION_STRINGENCY LENIENT
		fi
	elif [ $repeat_mapping = "no" ]; then

		# in case that one wants to repeat certain steps with the existing BAM files, we do not want to copy
		if [ $runmode = 'MS' ]; then

			if [ ! -f "$name/results/bam/$name.Normal.bam" ]; then
				cp $bam_normal $name/results/bam/$name.Normal.bam
				samtools index -@ $threads $name/results/bam/$name.Normal.bam
				cp $name/results/bam/$name.Normal.bai $name/results/bam/$name.Normal.bam.bai
			fi

			if [ ! -f "$name/results/bam/$name.Tumor.bam" ]; then
				cp $bam_tumor $name/results/bam/$name.Tumor.bam
				samtools index -@ $threads $name/results/bam/$name.Tumor.bam
				cp $name/results/bam/$name.Tumor.bai $name/results/bam/$name.Tumor.bam.bai
			fi

		elif [ $runmode = 'SS' ] && [ $types = 'Tumor' ]; then

			if [ ! -f "$name/results/bam/$name.Tumor.bam" ]; then
				cp $bam_tumor $name/results/bam/$name.Tumor.bam
				samtools index -@ $threads $name/results/bam/$name.Tumor.bam
				cp $name/results/bam/$name.Tumor.bai $name/results/bam/$name.Tumor.bam.bai
			fi

		elif [ $runmode = 'SS' ] && [ $types = 'Normal' ]; then

			if [ ! -f "$name/results/bam/$name.Normal.bam" ]; then
				cp $bam_normal $name/results/bam/$name.Normal.bam
				samtools index -@ $threads $name/results/bam/$name.Normal.bam
				cp $name/results/bam/$name.Normal.bai $name/results/bam/$name.Normal.bam.bai
			fi

		fi
	fi
fi


# 2. remap fastq
if [ $repeat_mapping = "yes" ]; then

	checkR1fq=$temp_dir/$name.${types}.R1.passed.fastq.gz
	checkR2fq=$temp_dir/$name.${types}.R2.passed.fastq.gz
	if [ -f "${checkR1fq}" ] && [ -f "${checkR2fq}" ]; then
		echo -e "skipping fastQC and trimming steps, found pre-existing files:\n ${checkR1fq}\n${checkR2fq}"
	else
		# echo '---- Calculating md5-sums ----' | tee -a $name/results/QC/$name.report.txt
		# echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
		# for type in $types;
		# do
		# md5sum $name/fastq/$name.$type.R1.fastq.gz > $name/fastq/$name.$type.R1.fastq.gz.md5 & PIDS="$PIDS $!"
		# md5sum $name/fastq/$name.$type.R2.fastq.gz > $name/fastq/$name.$type.R2.fastq.gz.md5 & PIDS="$PIDS $!"
		# done

		# wait $PIDS
		# PIDS=""

		echo '---- Running FastQC before trimming ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		for type in $types;
		do
		fastqc -t $threads \
		$name/fastq/$name.$type.R1.fastq.gz \
		$name/fastq/$name.$type.R2.fastq.gz \
		--outdir=$name/results/QC & PIDS="$PIDS $!"
		done

		wait $PIDS
		PIDS=""

		echo '---- Trimming reads ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		# FastQC has to be executed! all_DeterminePhred.sh will grep for the score in the QC output files
		for type in $types;
		do

		# check if the needed file exists
		fastqc_file=$name/results/QC/$name.$type.R1_fastqc.zip
		if [[ ! -f "$fastqc_file" ]]; then
		echo "FastQC file not found! Check if FastQC finished successfully and generated this file: $fastqc_file"
		exit 1
		fi

		trimmomatic_file=$(basename $trimmomatic_dir)
		if [ -z $phred ]; then phred=$(sh $repository_dir/all_DeterminePhred.sh $name $type); fi
		echo "Determined phred score for trimming: $phred"

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
		ILLUMINACLIP:$trimmomatic_dir/adapters/TruSeq3-PE-2.fa:2:30:10 & PIDS="$PIDS $!"
		done
		# output (in temp dir): _not_passed.fastq.gz, passed.fastq.gz for each type

		wait $PIDS
		PIDS=""

		echo '---- Running FastQC after trimming ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		for type in $types;
		do
		fastqc -t $threads \
		$temp_dir/$name.$type.R1.passed.fastq.gz \
		$temp_dir/$name.$type.R2.passed.fastq.gz \
		--outdir=$name/results/QC & PIDS="$PIDS $!"
		done

		wait $PIDS
		PIDS=""

		# echo '---- Removing fastq files ----' | tee -a $name/results/QC/$name.report.txt
		# echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		# moved to end of pipeline
		# find $name/fastq/ -type f -name "$name*fastq.gz" -exec rm -r {} + & PIDS="$PIDS $!"

		# wait $PIDS
		# PIDS=""
	fi
	
	
	
	checkSam=$temp_dir/$name.$type.mapped.sam
	if [ -f "${checkSam}" ]; then
		echo -e "skipping mapping step, found pre-existing file:\n ${checkSam}"
	else
	
		echo '---- Mapping trimmed reads ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
		
		type=$types

		bwa mem -t $threads $genomeindex_dir \
		-Y -K $bwainputbases -v 1 \
		$temp_dir/$name.$type.R1.passed.fastq.gz \
		$temp_dir/$name.$type.R2.passed.fastq.gz \
		-O $temp_dir/$name.$type.mapped.sam
	
	fi
	
	checkBam=$temp_dir/$name.${types}.cleaned.bam
	if [ -f "${checkBam}" ]; then
		echo -e "skipping mapping step, found pre-existing file:\n ${checkBam}"
	else
		echo '---- Cleaning mapped sam ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
		
		type=$types
		
		java -Xmx${RAM}G -Dpicard.useLegacyParser=false \
		-jar $picard_dir/picard.jar CleanSam \
		-I $temp_dir/$name.$type.mapped.sam \
		-O $temp_dir/$name.$type.cleaned.bam \
		-VALIDATION_STRINGENCY LENIENT
	
	fi
	
	checkBam=$temp_dir/$name.${types}.cleaned.sorted.readgroups.marked.bam
	if [ -f "${checkBam}" ]; then
		echo -e "skipping Postprocessing I (Sorting, fixing read groups and marking duplicates), found pre-existing files:\n ${checkBam}"
	else
		# remove all fastqs based on runname + fastq.gz (should be passed and not_passed from trimmomatic in between)
		# find $temp_dir -type f -name "$name*fastq.gz" -exec rm -r {} + # moved to end of pipeline

		echo '---- Postprocessing I (Sorting, fixing read groups and marking duplicates) ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		for type in $types;
		do
			echo "Postprocessing I: Marking duplicates for $type" | tee -a $name/results/QC/$name.report.txt &&

			# this will take the mapped bam, sort it by QNAME, convert to SAM
			# next give it to SAMBLASTER for read dup marking, followed by sorting back to coordinates and format to BAM
			/opt/bin/sambamba sort \
			--sort-by-name \
			-t $threads -m ${RAM}GB --tmpdir $temp_dir \
			-o /dev/stdout $temp_dir/$name.$type.cleaned.bam | samtools view -h | /opt/samblaster-0.1.26/samblaster | samtools view -Sb | /opt/bin/sambamba sort \
			-t $threads -m ${RAM}GB --tmpdir $temp_dir \
			-o $temp_dir/$name.$type.cleaned.sorted.marked.bam /dev/stdin &&

			# moved to end of pipeline
			# rm -f $temp_dir/$name.$type.cleaned.bam &&

			echo "Postprocessing I: Adding read groups for $type" | tee -a $name/results/QC/$name.report.txt &&
			java -Xmx${RAM}G -Dpicard.useLegacyParser=false \
			-jar $picard_dir/picard.jar AddOrReplaceReadGroups \
			-I $temp_dir/$name.$type.cleaned.sorted.marked.bam \
			-O $temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam \
			-ID 1 -LB Lib1 -PL ILLUMINA -PU Run1 -SM $type \
			-MAX_RECORDS_IN_RAM $MAX_RECORDS_IN_RAM & PIDS="$PIDS $!"

			# moved to end of pipeline
			# rm -f $temp_dir/$name.$type.cleaned.sorted.bam &&
			# rm -f $temp_dir/$name.$type.cleaned.sorted.bam.bai &&
			# rm -f $temp_dir/$name.$type.cleaned.sorted.marked.bam &&
			# rm -f $temp_dir/$name.$type.cleaned.sorted.marked.bam.bai &&
			# rm -f $temp_dir/$name.$type.cleaned.sorted.readgroups.bam & PIDS="$PIDS $!"
		done

		wait $PIDS
		PIDS=""

		# sambamba tends to break occasionally, so we check the files to stop Mutect2 from working on broken files
		for type in $types;
		do
			CheckFile $temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam
		done
	
	fi

	checkRecalTable=$name/results/QC/$name.${types}.GATK4.post.recal.table
	if [ -f "${checkRecalTable}" ]; then
		echo -e "skipping postprocessing II (Base recalibration), found pre-existing files:\n ${checkRecalTable}"
	else
	
		echo '---- Postprocessing II (Base recalibration) ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		for type in $types;
		do
		echo "Postprocessing II for $type"  | tee -a $name/results/QC/$name.report.txt &&
		java -Xmx${RAM}G -jar $GATK_dir/gatk.jar BaseRecalibrator \
		-R $genome_file \
		-I $temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam \
		--known-sites $snp_file \
		--use-original-qualities \
		-O $name/results/QC/$name.$type.GATK4.pre.recal.table &&

		# this will also generate a bam index (.bai)
		java -Xmx${RAM}G -jar $GATK_dir/gatk.jar ApplyBQSR \
		-R $genome_file \
		-I $temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam \
		-O $name/results/bam/$name.$type.bam \
		-bqsr $name/results/QC/$name.$type.GATK4.pre.recal.table &&

		# moved to end of pipeline
		# rm -f $temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam &&
		# rm -f $temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam.bai &&

		java -Xmx${RAM}G -jar $GATK_dir/gatk.jar BaseRecalibrator \
		-R $genome_file \
		-I $name/results/bam/$name.$type.bam \
		--known-sites $snp_file \
		--use-original-qualities \
		-O $name/results/QC/$name.$type.GATK4.post.recal.table &&

		cp $name/results/bam/$name.$type.bai $name/results/bam/$name.$type.bam.bai & PIDS="$PIDS $!"
		done

		wait $PIDS
		PIDS=""
	
	fi
fi

rm -rf '?'

echo '---- removing intermediate files ----'
find $name/fastq/ -type f -name "$name*fastq.gz" -exec rm -r {} +
find $temp_dir -type f -name "$name*fastq.gz" -exec rm -r {} +
# rm -f $temp_dir/$name.$types.cleaned.bam
# rm -f $temp_dir/$name.$type.mapped.sam
# rm -f $temp_dir/$name.$types.cleaned.sorted.bam
# rm -f $temp_dir/$name.$types.cleaned.sorted.bam.bai
# rm -f $temp_dir/$name.$types.cleaned.sorted.marked.bam
# rm -f $temp_dir/$name.$types.cleaned.sorted.marked.bam.bam
# rm -f $temp_dir/$name.$types.cleaned.sorted.readgroups.bam
# rm -f $temp_dir/$name.$types.cleaned.sorted.readgroups.marked.bam
# rm -f $temp_dir/$name.$types.cleaned.sorted.readgroups.marked.bam.bai

echo '---- Finished analysis of sample '$name' ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

exit 0
