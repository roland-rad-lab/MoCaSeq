# FORMAT INPUT FILES
# setwd("/mnt/3TBVol/MoCaSeq_WD/MoCaSeq_ref/GRCh38.p12/vcf_annotation/")
# dt <- fread("hg38_ucsc_repeatMasker.txt")
# dt <- dt[, .(genoName, genoStart, genoEnd, repName, swScore, strand)]
# dt[, genoName := gsub("chr", "", genoName)]
# dt <- dt[!grep("Un|random|GL|JH|M|KQ|KZ|KI|EBV|Phix|KN|KV|alt", genoName)]
# fwrite(dt, "hg38_repeatMasker.bed", sep="\t", col.names = F)

# agilent probes
# dt[, V1 := gsub("chr", "", V1)]
# fwrite(dt, "GRCh38_SureSelect.S31285117_MergedProbes_All_Exons_V7.bed", sep="\t", col.names = F)


# RUN THIS WITHIN DOCKER
# sudo docker run \
# -it --entrypoint=/bin/bash \
# -v ${working_directory}:/var/pipeline/ \
# -v ${ref_directory}:/var/pipeline/ref/ \
# -v ${script_directory}:/opt/MoCaSeq \
# -v ${working_directory}/temp/:/var/pipeline/temp/ \
# mocaseq2

# USING THIS COMMAND
#bash /opt/MoCaSeq/repository/vcf_annotation/Annotate_Variants_Wrapper.sh hPDAC03_HD_LivMet-1 Human

name=$1
species=$2
echo "RUNNING $name for species $species"



temp_dir=/var/pipeline/temp2
annTmp=${temp_dir}/annInput
mkdir -p ${annTmp}
# example: /home/rad/users/gaurav/projects/pdacMetastasis/output/mPDAC_wgs/annInput


genomes_dir=/var/pipeline/ref

# set some species specific arguments
if [ $species = 'Mouse' ]; then
	echo 'Species set to Mouse'
  genome_dir=$genomes_dir/GRCm38.p6
  genome_file=$genome_dir/GRCm38.p6.fna
  segdupbed=$genome_dir"/vcf_annotation/mm10_ucsc_segmentalDups.bed"
  repmaskbed=$genome_dir"/vcf_annotation/mm10_repeatMasker.bed"
  gcbed=$genome_dir"/vcf_annotation/mm10_5bp_GC_percentage.bedgraph"
  probesbed=$genome_dir"/GRCm38.SureSelect_Mouse_All_Exon_V1.bed"
elif [ $species = 'Human' ]; then
	echo 'Species set to Human'
  genome_dir=$genomes_dir/GRCh38.p12
  genome_file=$genome_dir/GRCh38.p12.fna
  segdupbed=$genome_dir"/vcf_annotation/hg38_ucsc_segmentalDups.bed"
  repmaskbed=$genome_dir"/vcf_annotation/hg38_repeatMasker.bed"
  gcbed=$genome_dir"/vcf_annotation/hg38_5bp_GC_percentage.bedgraph"
  probesbed=$genome_dir"/GRCh38_SureSelect.S31285117_MergedProbes_All_Exons_V7.bed"
else echo "Invalid species input (${species}). Choose Mouse or Human"; exit 1
fi

#genomeFasta="${jobdir}/input/annotation/${species}/GRCm38.p6.fna" # given by genome_file

# Input parameters
projname="${name}_annPipeline" #'mPDAC_wgs' # Just some name... you can use 'tmp'
#jobdir='/home/rad/users/gaurav/projects/pdacMetastasis' # do we even need this?!

annoScriptDir=/opt/MoCaSeq/repository/vcf_annotation/
# THIS IS ALSO HARDCODED WITHIN annotate_variant_bash_wrapper.py

tumbamFile=${name}/results/bam/${name}.Tumor.bam
norbamFile=${name}/results/bam/${name}.Normal.bam


# 0) Merge rescued files and save to temp
Rscript ${annoScriptDir}R_merge_vcf_rescued_files.R ${name} ${annTmp}
mutationfile=${annTmp}/${name}_Mutect2.txt


##ONLY NEEDED FOR NON RESCUE FILES
#mutationfile=${name}/results/Mutect2/${name}.Mutect2.txt
# example: /media/rad/HDD1/pdacMetastasis/input/mPDAC_wgs/DS4/results/Mutect2/DS4.Mutect2.txt
# 1) Convert the mutect2 format to pipeline mutect2 format
# Input:  Mutect2.txt
# Output: Mutect2.txt with additional columns in a new output folder
#python3.7 ${annoScriptDir}convert_vcf_to_annotation_input.py -if=${mutationfile} -od=${annTmp} -sn=${name}


# 2) Generate the bash wrapper
# This will generate a final bash wrapper in the script folder that is created in the annotate_variant_bash_wrapper.py file
python3.7 ${annoScriptDir}annotate_variant_bash_wrapper.py -gf=${genome_file} -sd=${segdupbed} -rp=${repmaskbed} -gc=${gcbed} -pb=${probesbed} -tb=${tumbamFile} -nb=${norbamFile} -pn=${projname} -pd=${annTmp} -sn=${name}

# Run the wrapper
mkdir -p ${annTmp}/finalAnnotation/
bash ${annTmp}/scripts/wrapper/${projname}/${name}_Mutect2_annotation_wrapper.sh

mkdir ${name}/results/VariantAnnotation
chown -R 1000:1000 ${name}/results/VariantAnnotation

# copy result files
echo "Copy result files"
cp ${annTmp}/finalAnnotation/*  ${name}/results/VariantAnnotation/
chown -R 1000:1000 ${name}/results/VariantAnnotation/${name}_Mutect2_variant_annotation.*

# copy log file
echo "Copy log files"
cp ${annTmp}/logs/${name}_annotation.log  ${name}/results/VariantAnnotation/
chown -R 1000:1000 ${name}/results/VariantAnnotation/${name}_annotation.log

# copy processed rescue file
echo "Copy combined rescue file"
cp ${annTmp}/${name}_Mutect2.txt  ${name}/results/rescued/${name}_Mutect2_combined.tsv
chown -R 1000:1000 ${name}/results/rescued/${name}_Mutect2_combined.tsv

rm -r ${annTmp} ${temp_dir}/interimFiles
echo "FINISHED"
