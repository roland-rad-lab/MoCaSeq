#!/bin/bash

##########################################################################################
##
## SNV_Mutect2Postprocessing.sh
##
## Postprocessing for Mutect2.
##
##########################################################################################

name=$1
runmode=$2
sequencing_type=$3
config_file=$4
species=$5
threads=$6
types=$7

. $config_file

AccessBed=${temp_dir}/${name}_access-CNVKit.bed
CNVKit_folder=$name/results/CNVKit
bam_tumor=$name/results/bam/$name.Tumor.bam
bam_normal=$name/results/bam/$name.Normal.bam

log_file=$name/results/QC/$name.report.txt
#log_file=$name/results/CNVKit/$name.debug.report.txt

# Run CNVKit for multisample for WES or WGS
# the commands for WGS/WES are identical, except for: -method, -access, -targets
if [ $runmode = "MS" ] && [ $sequencing_type = 'WES' ]; then
  echo '---- Running CNVKit (matched tumor-normal, WES) ----' | tee -a ${log_file}
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a ${log_file}

	mkdir ${CNVKit_folder}/matched/

  cnvkit.py access $genome_file -o ${AccessBed}

  cnvkit.py batch \
    $bam_tumor \
    --normal $bam_normal \
    --fasta "$genome_file" \
    --output-reference ${CNVKit_folder}/matched/Refernce.cnn \
    --output-dir ${CNVKit_folder}/matched/  \
    --short-names \
    --diagram \
    --scatter \
    --annotate "$RefFlat" \
    --access "$AccessBed" \
    --targets "$exons_file" \
    --drop-low-coverage \
    -m hybrid \
    -p "$threads"

    # maybe improve this some day:
    #  --antitarget-avg-size 50000 --> CNVkit uses a cautious default off-target bin size that, in our experience, will typically include more reads than the average on-target bin. However, we encourage the user to examine the coverage statistics reported by CNVkit and specify a properly calculated off-target bin size for their samples in order to maximize copy number information.

elif [ $runmode = "MS" ] && [ $sequencing_type = 'WGS' ]; then

  echo '---- Running CNVKit (matched tumor-normal, WGS) ----' | tee -a ${log_file}
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a ${log_file}

	mkdir ${CNVKit_folder}/matched/

  #To speed up and/or improve the accuracy of WGS analyses: instead of analyzing the whole genome, use the “target” BED file to limit the analysis to just the genic regions. (like described in the CNVKit vignette)
  cnvkit.py batch \
    $bam_tumor \
    --normal $bam_normal \
    --fasta "$genome_file" \
    --output-reference ${CNVKit_folder}/matched/Refernce.cnn \
    --output-dir ${CNVKit_folder}/matched/ \
    --short-names \
    --diagram \
    --scatter \
    --annotate "$RefFlat" \
    --access "" \
    --targets "$genecode_file_genes_bed" \
    --drop-low-coverage \
    -m wgs \
    -p "$threads"
fi

# rename the matched files to prevent any confusion
# rename the matched files to prevent any confusion
if [ $runmode = "MS" ]; then
cnvkit_files=$(find ${CNVKit_folder}/matched/ -type f -name "*Tumor*")

for file in $cnvkit_files
do
	mv "$file" "`echo $file | sed 's/Tumor.//'`"
done

fi


# Run CNVKit for single sample (and always after MS) for WES or WGS
echo '---- Running CNVKit (single-sample) ----' | tee -a ${log_file}

mkdir ${CNVKit_folder}/single/

#for single sample the normal is "empty"
for type in $types;
do

	echo "---- Running CNVKit (${type}) ----" | tee -a ${log_file}
	echo -e "$(date) \t timestamp: $(date +%s)" | tee -a ${log_file}

  if [ $sequencing_type = 'WES' ]; then
    cnvkit.py access $genome_file -o ${AccessBed}

    cnvkit.py batch \
      $name/results/bam/$name.$type.bam \
      --normal \
      --fasta "$genome_file" \
      --output-reference ${CNVKit_folder}/single/Refernce.${type}.cnn \
      --output-dir ${CNVKit_folder}/single/ \
      --short-names \
      --diagram \
      --scatter \
      --annotate "$RefFlat" \
      --access "$AccessBed" \
      --targets "$exons_file" \
      --drop-low-coverage \
      -m hybrid \
      -p "$threads"

  elif [ $sequencing_type = 'WGS' ]; then

    cnvkit.py batch \
      $name/results/bam/$name.$type.bam \
      --normal \
      --fasta "$genome_file" \
      --output-reference ${CNVKit_folder}/single/Refernce.${type}.cnn \
      --output-dir ${CNVKit_folder}/single/ \
      --short-names \
      --diagram \
      --scatter \
      --annotate "$RefFlat" \
      --access "" \
      --targets "$genecode_file_genes_bed" \
      --drop-low-coverage \
      -m wgs \
      -p "$threads"
  fi
done
