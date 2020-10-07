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
	echo "	-art, --artefact         Set to 'GT' (oxidation artefact), 'CT' (FFPE artefact) or 'none'. Optional. If set to something other than 'none' AND Mutect2 is 'yes', forces quality_control to 'yes'. Defaults to none."
	echo "	-filt, --filtering       Set to 'all' (AF >= 0.05, , Variant in Tumor >= 2, Variant in Normal <= 1, Coverage >= 5), 'hard' (AF >= 0.1, Variant in Tumor >= 3, Variant in Normal = 0, Coverage >= 10) or 'none' (no filters). Optional. Defaults to 'hard'."
	echo "	-p, --phred              If not set, script will try to automatically extract phred-score. Otherwise, set manually to 'phred33' or 'phred64'. 'phred64' only relevant for Illumina data originating before 2011. Optional."
	echo "	-mu, --Mutect2           Set to 'yes' or 'no'. Needed for LOH analysis and Titan. Greatly increases runtime for WGS. Optional. Defaults to 'yes'."
	echo "	-de, --Delly             Set to 'yes' or 'no'. Needed for chromothripsis inference. Do not use for WES. Optional. Defaults to 'no'. Only use in matched-sample mode."
	echo "	-ti, --Titan             Set to 'yes' or 'no'. Runs TITAN to model subclonal copy number alterations, predict LOH and estimate tumor purity. Greatly increases runtime for WGS. If set to 'yes', forces Mutect2 to 'yes'. Optional. Defaults to 'yes' for WES and 'no' for WGS. Only use in matched-sample mode."
	echo "	-abs, --Absolute         Set to 'yes' or 'no'. Runs ABSOLUTE to estimate purity/ploidy and compute copy-numbers. Optional. Can also include information from somatic mutation data, for this set Mutect2 to 'yes'."
	echo "	-fac, --Facets           Set to 'yes' or 'no'. Runs the allele-specific copy-number caller FACETS with sample purity estimations. Optional. Defaults to 'yes' for WES and 'no' for WGS. Only use in matched-sample mode."
	echo "	-bt, --BubbleTree             Set to 'yes' or 'no'. Runs the analysis of tumoral aneuploidy and clonality. Optional. If set to 'yes', forces Mutect2 to 'yes'. Optional. Defaults to 'yes' for WES and 'no' for WGS. Only use in matched-sample mode."
	echo "	-gatk, --GATKVersion     Set to '4.1.0.0', '4.1.3.0', '4.1.4.1' or '4.1.7.0', determining which GATK version is used. Optional. Defaults to 4.1.7.0"
	echo "	--test                   If set to 'yes': Will download reference files (if needed) and start a test run. All other parameters will be ignored"
	echo "	--memstats               If integer > 0 specified, will write timestamped memory usage and cumulative CPU time usage of the docker container to ./results/memstats.txt every <integer> seconds. Defaults to '0'."
	echo "	--help                   Show this help."
  exit 1
}

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
artefact_type=none
filtering=hard
phred=
Mutect2=yes
Titan=no
Absolute=no
Facets=no
BubbleTree=no
Delly=no
runmode=MS
GATK=4.1.7.0
test=no
memstats=0
config_file=
species=Mouse

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
    -de|--Delly) shift;Delly="$1";;
    -ti|--Titan) shift;Titan="$1";;
		-abs|--Absolute) shift;Absolute="$1";;
		-fac|--Facets) shift;Facets="$1";;
		-bt|--BubbleTree) shift;BubbleTree="$1";;
    -gatk|--GATKVersion) shift;GATK="$1";;
    --memstats) shift;memstats="$1";;
    --test) shift;test="$1";;
    --help) usage;shift;;
	*) usage;shift;;
esac; shift; done

if [ -z $config_file ]; then
	config_file=/opt/MoCaSeq/config.sh
fi

# init variables of test run
test_dir=${config_file%/*}/test

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
	artefact_type=none
	filtering=all
	phred=
	Mutect2=yes
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

# CHECK IF ALL NEEDED FILE EXIST
if [ "$fastq_normal_1" != "" ] && [ ! -f "$fastq_normal_1" ]; then
echo "ERROR, File not found: $fastq_normal_1"
exit 1
fi

if [ "$fastq_normal_2" != "" ] && [ ! -f "$fastq_normal_2" ]; then
echo "ERROR, File not found: $fastq_normal_2"
exit 1
fi

if [ "$fastq_tumor_1" != "" ] && [ ! -f "$fastq_tumor_1" ]; then
echo "ERROR, File not found: $fastq_tumor_1"
exit 1
fi

if [ "$fastq_tumor_2" != "" ] && [ ! -f "$fastq_tumor_2" ]; then
echo "ERROR, File not found: $fastq_tumor_2"
exit 1
fi

if [ "$bam_normal" != "" ] && [ ! -f "$bam_normal" ]; then
echo "ERROR, File not found: $bam_normal"
exit 1
fi

if [ "$bam_tumor" != "" ] && [ ! -f "$bam_tumor" ]; then
echo "ERROR, File not found: $bam_tumor"
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
	Titan=no

	if [ $Titan = 'yes' ]; then
	echo 'PARAMETER CHANGE: TITAN needs matched tumor/normal and can not be used on single sample runs. TITAN was set to "no".'
	fi
fi

# TODO: artifact will never be no since the default is artefact_type=none
if [ $Mutect2 = 'yes' ] && [ $artefact_type != 'none' ]; then
	quality_control=yes
	echo 'PARAMETER CHANGE: Quality control was set to "yes", due to mutect2="yes" and artifact!="none".'
else
	Titan=no

	if [ $Titan = 'yes' ]; then
	echo 'PARAMETER WARNING: TITAN needs Mutect2="yes" and artefact!="none". TITAN was set to "no".'
	fi
fi


# Absolute ist default 'yes' for WES and default 'no' for WGS
if [ $sequencing_type = 'WES' ] && [ -z $Absolute ]; then
	Absolute=yes
elif [ $sequencing_type = 'WGS' ] && [ -z $Absolute ]; then
	Absolute=no
fi
# Absolute needs segmented copy ratios data data ((/Copywriter/CNAprofiles/segment.Rdata or HMMCopy/HMMCopy.20000.segments.txt"))
# Absolute can optionally use somatic mutation data (Mutect2.vep.maf.fn), this will be determined by the value of Mutect2

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

# BubbleTree ist default 'yes' for WES and default 'no' for WGS
if [ $sequencing_type = 'WES' ] && [ -z $BubbleTree ]; then
	BubbleTree=yes
elif [ $sequencing_type = 'WGS' ] && [ -z $BubbleTree ]; then
	BubbleTree=no
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
	BubbleTree=no

	if [ $BubbleTree = 'yes' ]; then
	echo 'PARAMETER CHANGE: BubbleTree needs matched tumor/normal and can not be used on single sample runs. BubbleTree was set to "no".'
	fi
fi





#reading configuration from $config_file
source $config_file
repository_dir=${config_file%/*}/repository

echo "---- Starting ${species} Cancer Genome Analysis ----"
echo -e "$(date) \t timestamp: $(date +%s)"

echo '---- Creating directories ----'
echo -e "$(date) \t timestamp: $(date +%s)"
mkdir -p $name/
mkdir -p $name/pipeline/
mkdir -p $name/fastq/
mkdir -p $name/results/QC
mkdir -p $name/results/bam
mkdir -p $name/results/Manta
mkdir -p $name/results/Strelka
mkdir -p $name/results/msisensor

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
if [ ! grep -Fxq "DONE" $genome_dir/GetReferenceData.txt ]; then
	echo '---- Reference files not found - Files will be downloaded ----' | tee -a $name/results/QC/${name}.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/${name}.report.txt
  rm -rf $genome_dir
	if [ $species = 'Mouse' ]; then
	sh $repository_dir/Preparation_GetReferenceDataMouse.sh $config_file $temp_dir
	elif [ $species = 'Human' ]; then
	sh $repository_dir/Preparation_GetReferenceDataHuman.sh $config_file $temp_dir
	fi
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

if [ $sequencing_type = 'WES' ]; then
	mkdir -p $name/results/Copywriter
fi

if [ $runmode = 'MS' ]; then
	mkdir -p $name/results/HMMCopy
fi

if [ $Mutect2 = 'yes' ]; then
	mkdir -p $name/results/Mutect2
fi

if [ $Mutect2 = 'yes' ] && [ $runmode = 'MS' ]; then
	mkdir -p $name/results/LOH
fi

if [ $Titan = 'yes' ] && [ $Mutect2 = 'yes' ] && [ $runmode = 'MS' ]; then
	mkdir -p $name/results/Titan
fi

if [ $Absolute = 'yes' ]; then
	mkdir -p $name/results/ABSOLUTE
fi

if [ $Facets = 'yes' ]; then
	mkdir -p $name/results/FACETS
fi

if [ $BubbleTree = 'yes' ]; then
	mkdir -p $name/results/BubbleTree
fi

if [ $Delly = 'yes' ] && [ $runmode = 'MS' ] && [ $sequencing_type = 'WGS' ]; then
	mkdir -p $name/results/Delly
	mkdir -p $name/results/Chromothripsis
fi

if [ $species = 'Mouse' ]; then
	mkdir -p $name/results/Genotype
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
echo Assuming $artefact_type-artefacts for SNV-calling | tee -a $name/results/QC/$name.report.txt
echo $filtering is setting for filtering of SNV calls | tee -a $name/results/QC/$name.report.txt

echo Quality scores are assumed as $phred | tee -a $name/results/QC/$name.report.txt

echo Using GATK v$GATK | tee -a $name/results/QC/$name.report.txt
if [ $Mutect2 = "yes" ]; then
	echo Will run Mutect2 | tee -a $name/results/QC/$name.report.txt
fi
if [ $Delly = "yes" ]; then
	echo Will run Delly | tee -a $name/results/QC/$name.report.txt
fi
if [ $Titan = "yes" ]; then
	echo Will run Titan | tee -a $name/results/QC/$name.report.txt
fi
if [ $Absolute = "yes" ]; then
	echo Will run ABSOLUTE | tee -a $name/results/QC/$name.report.txt
fi
if [ $Facets = "yes" ]; then
	echo Will run FACETS | tee -a $name/results/QC/$name.report.txt
fi
if [ $BubbleTree = "yes" ]; then
	echo Will run BubbleTree | tee -a $name/results/QC/$name.report.txt
fi
echo Starting workflow using $threads CPU-threads and $RAM GB of RAM | tee -a $name/results/QC/$name.report.txt




#rerouting STDERR to report file
exec 2>> $name/results/QC/$name.report.txt
#exec 2>&1 | tee $name/results/QC/$name.report.txt





echo '---- Creating directories ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

echo '---- Copying repository ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

cp -r $repository_dir/ $name/pipeline/
cp $repository_dir/../MoCaSeq.sh $name/pipeline/
cp $config_file $name/pipeline/

echo '---- Copying raw data ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

# 1. this will either copy fastqs to location OR do sam2fastq
if [ -z $bam_normal ] && [ -z $bam_tumor ]; then

	if [ $runmode = 'MS' ]; then
	#rsync /home/tin/file2.txt /home/tin/file3.txt /home/tin/Documents
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
	if [ $runmode = 'MS' ]; then
	cp $bam_normal $name/results/bam/$name.Normal.bam
	samtools index -@ $threads $name/results/bam/$name.Normal.bam
	cp $bam_tumor $name/results/bam/$name.Tumor.bam
	samtools index -@ $threads $name/results/bam/$name.Tumor.bam

	elif [ $runmode = 'SS' ] && [ $types = 'Tumor' ]; then
	cp $bam_tumor $name/results/bam/$name.Tumor.bam
	samtools index -@ $threads $name/results/bam/$name.Tumor.bam

	elif [ $runmode = 'SS' ] && [ $types = 'Normal' ]; then
	cp $bam_normal $name/results/bam/$name.Normal.bam
	samtools index -@ $threads $name/results/bam/$name.Normal.bam
	fi
fi

# 2. remap fastq
if [ $repeat_mapping = "yes" ]; then
	echo '---- Calculating md5-sums ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
	for type in $types;
	do
	md5sum $name/fastq/$name.$type.R1.fastq.gz > $name/fastq/$name.$type.R1.fastq.gz.md5 & PIDS="$PIDS $!"
	md5sum $name/fastq/$name.$type.R2.fastq.gz > $name/fastq/$name.$type.R2.fastq.gz.md5 & PIDS="$PIDS $!"
	done

	wait $PIDS
	PIDS=""

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

	echo '---- Removing fastq files ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	#TODO: SET RM
	#find $name/fastq/ -type f -name "$name*fastq.gz" -exec rm -r {} + & PIDS="$PIDS $!"

	wait $PIDS
	PIDS=""

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
	# output (in temp dir): .cleaned.bam for each type

	# remove all fastqs based on runname + fastq.gz (should be passed and not_passed from trimmomatic in between)
	#TODO: SET RM
	#find $temp_dir -type f -name "$name*fastq.gz" -exec rm -r {} + & PIDS="$PIDS $!"
	wait $PIDS
	PIDS=""

	echo '---- Postprocessing I (Sorting, fixing read groups and marking duplicates) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	for type in $types;
	do
	/opt/bin/sambamba sort \
	-t $threads -m ${RAM}GB --tmpdir $temp_dir \
	-o $temp_dir/$name.$type.cleaned.sorted.bam \
	$temp_dir/$name.$type.cleaned.bam #&&
	#NIKLAS TODO: sambamba can break for large files
	#TODO: SET RM

	#rm $temp_dir/$name.$type.cleaned.bam &&

	java -Xmx${RAM}G -Dpicard.useLegacyParser=false \
	-jar $picard_dir/picard.jar AddOrReplaceReadGroups \
	-I $temp_dir/$name.$type.cleaned.sorted.bam \
	-O $temp_dir/$name.$type.cleaned.sorted.readgroups.bam \
	-ID 1 -LB Lib1 -PL ILLUMINA -PU Run1 -SM $type \
	-MAX_RECORDS_IN_RAM $MAX_RECORDS_IN_RAM #&&

	#TODO: SET RM
	#rm $temp_dir/$name.$type.cleaned.sorted.bam &&
	#rm $temp_dir/$name.$type.cleaned.sorted.bam.bai &&

	/opt/bin/sambamba markdup --tmpdir $temp_dir \
	--t $threads \
	--overflow-list-size=$HASH_TABLE_SIZE --hash-table-size=$HASH_TABLE_SIZE \
	$temp_dir/$name.$type.cleaned.sorted.readgroups.bam \
	$temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam #&&

	#TODO: SET RM
	#rm $temp_dir/$name.$type.cleaned.sorted.readgroups.bam & PIDS="$PIDS $!"
	done

	wait $PIDS
	PIDS=""

	echo '---- Postprocessing II (Base recalibration) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt


	for type in $types;
	do
	echo "Running BaseRecalibrator 1 for $type"
	java -Xmx${RAM}G -jar $GATK_dir/gatk.jar BaseRecalibrator \
	-R $genome_file \
	-I $temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam \
	--known-sites $snp_file \
	--use-original-qualities \
	-O $name/results/QC/$name.$type.GATK4.pre.recal.table &&

	echo "Running ApplyBQSR for $type"
	java -Xmx${RAM}G -jar $GATK_dir/gatk.jar ApplyBQSR \
	-R $genome_file \
	-I $temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam \
	-O $name/results/bam/$name.$type.bam \
	-bqsr $name/results/QC/$name.$type.GATK4.pre.recal.table #&&

	#rm $temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam &&
	#rm $temp_dir/$name.$type.cleaned.sorted.readgroups.marked.bam.bai &&

	echo "Running BaseRecalibrator 2 for $type"
	java -Xmx${RAM}G -jar $GATK_dir/gatk.jar BaseRecalibrator \
	-R $genome_file \
	-I $name/results/bam/$name.$type.bam \
	--known-sites $snp_file \
	--use-original-qualities \
	-O $name/results/QC/$name.$type.GATK4.post.recal.table &&

	echo "Running Sambamba Index for $type"
	/opt/bin/sambamba index \
	-t $threads \
	$name/results/bam/$name.$type.bam \
	$name/results/bam/$name.$type.bam.bai #&&

	#TODO: SET RM
	#rm $name/results/bam/$name.$type.bai & PIDS="$PIDS $!"
	done

	wait $PIDS
	PIDS=""
fi

# RUN QUALITY CONTROL
if [ $quality_control = "yes" ]; then
	echo '---- Quality control I (Sequencing artifacts, multiple metrics) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	for type in $types;
	do
	java -Xmx${RAM}G -Dpicard.useLegacyParser=false \
	-jar $picard_dir/picard.jar CollectSequencingArtifactMetrics \
	-R $genome_file \
	-I $name/results/bam/$name.$type.bam \
	-O $name/results/QC/$name.$type.bam.artifacts & PIDS="$PIDS $!"

	java -Xmx${RAM}G -Dpicard.useLegacyParser=false \
	-jar $picard_dir/picard.jar CollectMultipleMetrics \
	-R $genome_file \
	-I $name/results/bam/$name.$type.bam \
	-O $name/results/QC/$name.$type.bam.metrics & PIDS="$PIDS $!"

	samtools idxstats $name/results/bam/$name.$type.bam \
	> $name/results/QC/$name.$type.bam.idxstats & PIDS="$PIDS $!"
	done

	wait $PIDS
	PIDS=""

	echo '---- Quality control II (WES- or WGS-specific metrics) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	if [ $sequencing_type = 'WES' ]; then
	for type in $types;
	do
		java -Xmx${RAM}G -Dpicard.useLegacyParser=false \
		-jar $picard_dir/picard.jar CollectHsMetrics \
		-SAMPLE_SIZE 100000 \
		-R $genome_file \
		-I $name/results/bam/$name.$type.bam \
		-O $name/results/QC/$name.$type.bam.metrics \
		-BAIT_INTERVALS $interval_file \
		-TARGET_INTERVALS $interval_file & PIDS="$PIDS $!"
	done

	wait $PIDS
	PIDS=""

	elif [ $sequencing_type = 'WGS' ]; then
	for type in $types;
	do
		java -Xmx${RAM}G -Dpicard.useLegacyParser=false \
		-jar $picard_dir/picard.jar CollectWgsMetrics \
		-R $genome_file \
		-I $name/results/bam/$name.$type.bam \
		-O $name/results/QC/$name.$type.bam.metrics \
		-SAMPLE_SIZE 1000000 & PIDS="$PIDS $!"
	done

	wait $PIDS
	PIDS=""
	fi

	echo '---- Summarizing quality control data ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	multiqc $name/results/QC -n $name -o $name/results/QC/ --pdf --interactive

fi

# CHECK FOR CORRECT BAM PAIRING
if [ $runmode = "MS" ]; then
	echo '---- Check for matched BAM-files ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	python2 $bammatcher_dir/bam-matcher.py \
	-B1 $name/results/bam/$name.Tumor.bam \
	-B2 $name/results/bam/$name.Normal.bam \
	--config $bammatcher_file --html --number_of_snps 100000 \
	--output $name/results/QC/$name.Tumor.Normal.bammatcher.txt
fi

# SET GENOTYPES FOR MOUSE
if [ $species = "Mouse" ]; then
	echo '---- Get genotypes ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	echo -e 'Name\tAllele\tCHROM\tPOS\tREF\tALT\tCount_Ref\tCount_Alt\tComment' > $name/results/Genotype/$name.Genotypes.txt

	allele=Kras-G12D
	position=6:145246771-145246771
	comment="GGT>GAT=G>D"
	sh $repository_dir/SNV_GetGenotype.sh $name $allele $comment $config_file $species $position $runmode $types

	allele=Kras-G12D_Neighbour
	position=6:145246791-145246791
	comment="shift_always_identical_to_Kras-G12D_20bp_upstream"
	sh $repository_dir/SNV_GetGenotype.sh $name $allele $comment $config_file $species $position $runmode $types

	allele=Trp53-R172H
	position=11:69588512-69588512
	comment="CGC>CAC=R>H"
	sh $repository_dir/SNV_GetGenotype.sh $name $allele $comment $config_file $species $position $runmode $types

	allele=BrafLSL-637E
	position=6:39627783-39627783
	comment="CAC>CTC=V>E"
	sh $repository_dir/SNV_GetGenotype.sh $name $allele $comment $config_file $species $position $runmode $types

	allele=Ink4aKO
	position=4:89276809-89276811
	comment="GCC>TAG=A>*"
	sh $repository_dir/SNV_GetGenotype.sh $name $allele $comment $config_file $species $position $runmode $types

	if [ $runmode = "MS" ]; then
	allele=Trp53-fl
	position=11:69580359-69591872
	transcript=ENSMUST00000171247.7
	wt_allele=1,11
	del_allele=2,3,4,5,6,7,8,9,10
	sh $repository_dir/CNV_GetGenotype.sh $name $position
	Rscript $repository_dir/CNV_GetGenotype.R $name $genecode_file_exons $transcript $allele $position $wt_allele $del_allele
	cat $name/results/Genotype/$name.Genotypes.temp.CNV.txt >> $name/results/Genotype/$name.Genotypes.txt
	rm $name/results/Genotype/$name.Genotypes.temp.CNV.txt

	allele=Cdh1-fl
	position=8:106603351-106670246
	transcript=ENSMUST00000000312.11
	wt_allele="1,2,3,16"
	del_allele="4,5,6,7,8,9,10,11,12,13,14,15"
	sh $repository_dir/CNV_GetGenotype.sh $name $position
	Rscript $repository_dir/CNV_GetGenotype.R $name $genecode_file_exons $transcript $allele $position $wt_allele $del_allele
	cat $name/results/Genotype/$name.Genotypes.temp.CNV.txt >> $name/results/Genotype/$name.Genotypes.txt
	rm $name/results/Genotype/$name.Genotypes.temp.CNV.txt

	allele=Pdk1-fl
	position=2:71873224-71903858
	transcript=ENSMUST00000006669.5
	wt_allele="1,2,5,6,7,8,9,10,11"
	del_allele="3,4"
	sh $repository_dir/CNV_GetGenotype.sh $name $position
	Rscript $repository_dir/CNV_GetGenotype.R $name $genecode_file_exons $transcript $allele $position $wt_allele $del_allele
	cat $name/results/Genotype/$name.Genotypes.temp.CNV.txt >> $name/results/Genotype/$name.Genotypes.txt
	rm $name/results/Genotype/$name.Genotypes.temp.CNV.txt

	allele=Raf1-fl
	position=6:115618569-115676635
	transcript=ENSMUST00000000451.13
	wt_allele="1,2,4,5,6,7,8,9,10,11,12,13,14,15,16,17"
	del_allele="3"
	sh $repository_dir/CNV_GetGenotype.sh $name $position
	Rscript $repository_dir/CNV_GetGenotype.R $name $genecode_file_exons $transcript $allele $position $wt_allele $del_allele
	cat $name/results/Genotype/$name.Genotypes.temp.CNV.txt >> $name/results/Genotype/$name.Genotypes.txt
	rm $name/results/Genotype/$name.Genotypes.temp.CNV.txt
	fi

	for type in $types;
	do
	rm $name/results/Genotype/$name.$type.Genotypes.CNV.txt
	done
	rm Rplots.pdf
fi


echo '---- Running Manta (matched tumor-normal) ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

if [ $sequencing_type = 'WES' ] && [ $runmode = "MS" ]; then
	python2 $manta_dir/bin/configManta.py \
	--normalBam $name/results/bam/$name.Normal.bam \
	--tumorBam $name/results/bam/$name.Tumor.bam \
	--referenceFasta $genome_file --runDir $name/results/Manta \
	--callRegions $callregions_file --generateEvidenceBam --exome
elif [ $sequencing_type = 'WGS' ] && [ $runmode = "MS" ]; then
	python2 $manta_dir/bin/configManta.py \
	--normalBam $name/results/bam/$name.Normal.bam \
	--tumorBam $name/results/bam/$name.Tumor.bam \
	--referenceFasta $genome_file --runDir $name/results/Manta \
	--callRegions $callregions_file --generateEvidenceBam
elif [ $sequencing_type = 'WES' ] && [ $runmode = "SS" ]; then
	python2 $manta_dir/bin/configManta.py \
	--bam $name/results/bam/$name.$types.bam \
	--referenceFasta $genome_file --runDir $name/results/Manta \
	--callRegions $callregions_file --generateEvidenceBam --exome
elif [ $sequencing_type = 'WGS' ] && [ $runmode = "SS" ]; then
	python2 $manta_dir/bin/configManta.py \
	--bam $name/results/bam/$name.$types.bam \
	--referenceFasta $genome_file --runDir $name/results/Manta \
	--callRegions $callregions_file --generateEvidenceBam
fi

python2 $name/results/Manta/runWorkflow.py -m local -j $threads

sh $repository_dir/SV_MantaPostprocessing.sh $name $species $config_file $runmode $types

if [ $runmode = "MS" ]; then
	sh $repository_dir/SNV_RunVEP.sh $name $config_file $species Manta MS
fi

# RUN STRELKA
echo '---- Running Strelka (matched tumor-normal) ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

if [ $sequencing_type = 'WES' ] && [ $runmode = "MS" ]; then
	python2 $strelka_dir/bin/configureStrelkaSomaticWorkflow.py \
	--normalBam $name/results/bam/$name.Normal.bam \
	--tumorBam $name/results/bam/$name.Tumor.bam \
	--ref $genome_file --runDir $name/results/Strelka/Strelka \
	--indelCandidates $name/results/Manta/results/variants/candidateSmallIndels.vcf.gz \
	--callRegions $callregions_file --exome
elif [ $sequencing_type = 'WGS' ] && [ $runmode = "MS" ]; then
	python2 $strelka_dir/bin/configureStrelkaSomaticWorkflow.py \
	--normalBam $name/results/bam/$name.Normal.bam \
	--tumorBam $name/results/bam/$name.Tumor.bam \
	--ref $genome_file --runDir $name/results/Strelka/Strelka \
	--indelCandidates $name/results/Manta/results/variants/candidateSmallIndels.vcf.gz \
	--callRegions $callregions_file
fi

if [ $runmode = "MS" ]; then
	python2 $name/results/Strelka/Strelka/runWorkflow.py -m local -j $threads

	sh $repository_dir/SNV_StrelkaPostprocessing.sh \
	$name $species $config_file $filtering $artefact_type $GATK

	Rscript $repository_dir/SNV_SelectOutput.R $name Strelka $species $CGC_file $TruSight_file

	sh $repository_dir/SNV_RunVEP.sh $name $config_file $species Strelka $runmode $types
fi

if [ $sequencing_type = 'WES' ]; then
	for type in $types;
	do
	python2 $strelka_dir/bin/configureStrelkaGermlineWorkflow.py \
	--bam $name/results/bam/$name.$type.bam \
	--ref $genome_file --runDir $name/results/Strelka/Strelka-$type \
	--callRegions $callregions_file --exome &&

	python2 $name/results/Strelka/Strelka-$type/runWorkflow.py -m local -j $threads & PIDS="$PIDS $!"
	done

	wait $PIDS
	PIDS=""
elif [ $sequencing_type = 'WGS' ]; then
	for type in $types;
	do
	python2 $strelka_dir/bin/configureStrelkaGermlineWorkflow.py \
	--bam $name/results/bam/$name.$type.bam \
	--ref $genome_file --runDir $name/results/Strelka/Strelka-$type \
	--callRegions $callregions_file &&

	python2 $name/results/Strelka/Strelka-$type/runWorkflow.py -m local -j $threads & PIDS="$PIDS $!"
	done

	wait $PIDS
	PIDS=""
fi

# RUN MUTECT2
# 1. matched tumor-normal
if [ $Mutect2 = 'yes' ] && [ $runmode = "MS" ]; then
	echo '---- Running Mutect2 (matched tumor-normal) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	java -Xmx${RAM}G -jar $GATK_dir/gatk.jar Mutect2 \
	--native-pair-hmm-threads $threads \
	-R $genome_file \
	-I $name/results/bam/$name.Tumor.bam \
	-I $name/results/bam/$name.Normal.bam \
	-tumor Tumor -normal Normal \
	--f1r2-tar-gz $name/results/Mutect2/$name.m2.f1r2.tar.gz \
	-O $name/results/Mutect2/$name.m2.vcf \
	-bamout $name/results/Mutect2/$name.m2.bam

	echo '---- Mutect2 Postprocessing (matched tumor-normal) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	sh $repository_dir/SNV_Mutect2Postprocessing.sh \
	$name $species $config_file $filtering $artefact_type $GATK

	Rscript $repository_dir/SNV_SelectOutput.R $name Mutect2 $species $CGC_file $TruSight_file

	Rscript $repository_dir/SNV_Signatures.R $name $species
fi

# 1. single-sample
if [ $Mutect2 = 'yes' ]; then
	echo '---- Running Mutect2 (single-sample) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	for type in $types;
	do
	java -Xmx${RAM}G -jar $GATK_dir/gatk.jar Mutect2 \
	--native-pair-hmm-threads $threads \
	-R $genome_file \
	-I $name/results/bam/$name.$type.bam \
	-tumor $type \
	--f1r2-tar-gz $name/results/Mutect2/$name."$type".m2.f1r2.tar.gz \
	-O $name/results/Mutect2/$name."$type".m2.vcf \
	-bamout $name/results/Mutect2/$name."$type".m2.bam

	sh $repository_dir/SNV_Mutect2PostprocessingSS.sh \
	$name $species $config_file $type $filtering $artefact_type $GATK

	Rscript $repository_dir/SNV_SelectOutputSS.R $name $type $species $CGC_file $TruSight_file
	done
fi

sh $repository_dir/SNV_RunVEP.sh $name $config_file $species Mutect2 $runmode $types

if [ $Mutect2 = 'yes' ] && [ $runmode = "MS" ]; then
	echo '---- Generate LOH data ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/LOH_GenerateVariantTable.R \
	$name $genome_file $repository_dir

	Rscript $repository_dir/LOH_MakePlots.R \
	$name $species $repository_dir
fi

echo '---- Generate and plot copy number data ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

resolution=20000

# RUN Copywriter
if [ $sequencing_type = 'WES' ]; then

	echo '---- Run CopywriteR ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/CNV_RunCopywriter.R $name $species $threads $runmode $genome_dir $centromere_file $varregions_file $resolution $types

	echo '---- Export raw data and re-normalize using Mode ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/CNV_CopywriterGetRawData.R $name $runmode $types
	python2 $repository_dir/CNV_CopywriterGetModeCorrectionFactor.py $name
	Rscript $repository_dir/CNV_CopywriterGetModeCorrectionFactor.R $name $runmode $types

	echo '---- Plot CNV-profiles ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/CNV_PlotCopywriter.R $name $species $repository_dir
	Rscript $repository_dir/CNV_MapSegmentsToGenes.R $name $species $genecode_file_genes Copywriter $resolution $CGC_file $TruSight_file

	# clean up
	find "$name/results/Copywriter/" -type f -name "$name*SegmentsChromosome*" -exec rm -r {} +
fi

# RUN HMMCopy (bin-size 20000)
if [ $runmode = "MS" ]; then
	echo '---- Run HMMCopy (bin-size 20000) ----' | tee -a $name/results/QC/$name.report.txt

	if [ $sequencing_type = 'WES' ]; then
	echo '     ---- HMMCopy is executed on this WES data only for testing purposes, please handle the results with caution! ----' | tee -a $name/results/QC/$name.report.txt
	fi

	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	sh $repository_dir/CNV_RunHMMCopy.sh $name $species $config_file $runmode $resolution $types
fi

# RUN HMMCopy (bin-size 1000)
if [ $runmode = "MS" ] && [ $sequencing_type = 'WGS' ]; then

	echo '---- Run HMMCopy (bin-size 1000) ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	sh $repository_dir/CNV_RunHMMCopy.sh $name $species $config_file $runmode 1000 $types
	Rscript $repository_dir/CNV_PlotHMMCopy.R $name $species $repository_dir $sequencing_type 1000 \
	$mapWig_file $gcWig_file $centromere_file $varregions_file
	Rscript $repository_dir/CNV_MapSegmentsToGenes.R $name $species $genecode_file_genes HMMCopy 1000 $CGC_file $TruSight_file
fi

# for MS + WES, PlotHMMCopy will be called but unable to find the HMMCopy/name.HMMCopy.20000.segments.txt which would be generated in the previous step (but is not because of sequencing_type=WGS)

# PLOT HMMCopy
if [ $runmode = "MS" ]; then

	echo '---- Plot HMMCopy ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/CNV_PlotHMMCopy.R $name $species $repository_dir $sequencing_type $resolution \
	$mapWig_file $gcWig_file $centromere_file $varregions_file
	Rscript $repository_dir/CNV_MapSegmentsToGenes.R $name $species $genecode_file_genes HMMCopy $resolution $CGC_file $TruSight_file
fi


echo '---- Run msisensor----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
if [ $runmode = "MS" ]; then
	msisensor msi -n $name/results/bam/$name.Normal.bam \
	-t $name/results/bam/$name.Tumor.bam \
	-o $name/results/msisensor/"$name".msisensor \
	-d $microsatellite_file -b $threads
elif [ $runmode = "SS" ]; then
	msisensor msi -t $name/results/bam/$name.$types.bam \
	-o $name/results/msisensor/"$name".$types.msisensor \
	-d $microsatellite_file -b $threads
fi

# PURITY ANALYSIS

# set CNV method (for ABSOLUTE and BubbleTree)
if [ $sequencing_type = 'WES' ]; then
CNVmethod="Copywriter"
elif [ $sequencing_type = 'WGS' ]; then
CNVmethod="HMMCopy"
fi

if [ $Titan = "yes" ]; then
	echo '---- Run TitanCNA ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/all_RunTitanCNA.R $name $species $repository_dir $resolution $mapWig_file $gcWig_file $exons_file $sequencing_type

	# rerun TITAN for each ploidy (2,3,4) and clusters (1 to numClusters)
	sh $repository_dir/all_RunTitanCNA.sh $name $repository_dir $threads $sequencing_type

	echo '		---- Mapping segments to genes ----' | tee -a $name/results/QC/$name.report.txt
	Rscript $repository_dir/LOH_MapSegmentsToGenes.R $name $species $genecode_file_genes $CGC_file $TruSight_file

	# cleanup TITAN temp dirs
	find . -maxdepth 1 -type d -name "run_ploidy*" -exec rm -r {} +
fi

if [ $Absolute = "yes" ]; then
	echo '---- Run ABSOLUTE ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	# Mutect2 can be "yes" or "no"
	Rscript $repository_dir/all_RunABSOLUTE.R $name $CNVmethod $species $Mutect2

	rm -r $name/results/ABSOLUTE/tmp/ # remove the tmp folder coming from R
fi


if [ $Facets = "yes" ]; then
	echo '---- Run FACETS ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	# first generate SNP pileup for FACETS
	echo '     ---- Run snp-pileup ----' | tee -a $name/results/QC/$name.report.txt
	bash ${repository_dir}/all_RunFACETS.sh $name $species $sequencing_type $config_file

	cval=300
	ndepth=10

	echo '     ---- Run FACETS analysis ----' | tee -a $name/results/QC/$name.report.txt
	Rscript $repository_dir/all_RunFACETS.R $name $species $cval $ndepth
fi

if [ $BubbleTree = "yes" ]; then
	echo '---- Run BubbleTree ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	Rscript $repository_dir/all_RunBubbleTree.R $name $CNVmethod
fi




if [ $sequencing_type = 'WGS' ] && [ $Delly = 'yes' ] && [ $runmode = "MS" ]; then
	echo '---- Optional for WGS: Run Delly ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	delly call \
	-o $name/results/Delly/$name.pre.bcf \
	-g $genome_file \
	$name/results/bam/$name.Tumor.bam \
	$name/results/bam/$name.Normal.bam

	delly filter \
	-f somatic -o $name/results/Delly/$name.bcf \
	-s $genome_dir/Samples.tsv $name/results/Delly/$name.pre.bcf

	bcftools view $name/results/Delly/$name.pre.bcf \
	> $name/results/Delly/$name.pre.vcf
fi

if [ $sequencing_type = 'WGS' ] && [ $Delly = 'yes' ] && [ $runmode = "MS" ]; then
	echo '---- Optional for WGS: Infer chromothripsis ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	format="tif"
	resolution=20000

	echo '---- Preparing input files and calculating coverage ----' | tee -a $name/results/QC/$name.report.txt
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

	coverage=$(sh $repository_dir/Chromothripsis_GetCoverage.sh $name)

	sh $repository_dir/Chromothripsis_FormatTable.sh $name

	Rscript $repository_dir/Chromothripsis_AnnotateRatios.R \
	-i $name/results/Delly/$name.breakpoints.tab \
	 > $name/results/Delly/$name.breakpoints_annotated.tab

	Rscript $repository_dir/Chromothripsis_FilterDelly.R \
	-n $name -c $coverage \
	-i $name/results/Delly/$name.breakpoints_annotated.tab

	for chr in $( seq $chromosomes ); do
	if [ $(Rscript $repository_dir/Chromothripsis_RearrangementCounter.R -i $name/results/Delly/$name.breakpoints.filtered.tab -c $chr) -ge 4 ]; then
		echo 'Analysing Chromosome '$chr
		mkdir -p $name'/results/Chromothripsis/Chr'$chr
		echo '---- Hallmark: Clustering of breakpoints for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_DetectBreakpointClustering.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-c $chr -n $name -f $format

		echo '---- Hallmark: Regularity of oscillating copy number states for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_SimulateCopyNumberStates.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-o $species_lowercase -c $chr -n $name -s 1000 -a 1000 -f $format -v 1

		echo '---- Hallmark: Interspersed loss and retention of heterozygosity for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_PlotLOHPattern.R \
		-s $name/results/HMMCopy/$name.HMMCopy.$resolution.segments.txt \
		-d $name/results/HMMCopy/$name.HMMCopy.$resolution.log2RR.txt \
		-v $name/results/LOH/$name.VariantsForLOH.txt \
		-o $species_lowercase -c $chr -n $name -f $format

		echo '---- Hallmark: Randomness of DNA fragment joins and segment order for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_DetectRandomJoins.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-c $chr -n $name -f $format

		echo '---- Hallmark: Ability to walk the derivative chromosome for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_WalkDerivativeChromosome.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-c $chr -n $name -f $format

		echo '---- Visualisation: Copy number profile combined with complex rearrangements for Chr'$chr' ----' | tee -a $name/results/QC/$name.report.txt
		echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

		Rscript $repository_dir/Chromothripsis_PlotRearrangementGraph.R \
		-i $name/results/Delly/$name.breakpoints.filtered.tab \
		-d $name/results/HMMCopy/$name.HMMCopy.$resolution.log2RR.txt \
		-c $chr -n $name -f $format
	else
		echo 'There are too few rearrangements in chromosome '$chr'.'
	fi
	done

fi

rm -rf '?'

echo '---- Finished analysis of sample '$name' ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

exit 0
