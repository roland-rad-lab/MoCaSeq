#!/bin/bash

##########################################################################################
##
## Preparation_GetReferenceData.sh
##
## Main routine for the download of all reference data needed for the WES and WGS workflows.
##
##########################################################################################

config_file=$1

#reading configuration from $config_file
. $config_file

VersionMouse=GRCm38.p6

mkdir Genomes
mkdir "Genomes/"$VersionMouse

echo '---- Get reference data ----' | tee "Genomes/"$VersionMouse"/GetReferenceData.txt"
echo Generate reference data for Version $VersionMouse | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"

#rerouting STDERR to report file
exec 2>> "Genomes/"$VersionMouse"/GetReferenceData.txt"

echo '---- Downloading reference genome ----' | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"

wget -P Genomes/"$VersionMouse" "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.8_"$VersionMouse"/GCA_000001635.8_"$VersionMouse"_genomic.fna.gz"
gunzip Genomes/"$VersionMouse"/"GCA_000001635.8_"$VersionMouse"_genomic.fna.gz"
mv Genomes/"$VersionMouse"/GCA_000001635.8_"$VersionMouse"_genomic.fna "Genomes/"$VersionMouse"/"$VersionMouse".fna"

echo '---- Generate BWA Index ----' | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"

sh $repository_dir/Preparation_GenerateBWAIndex.sh $VersionMouse $fasta_to_fastq

echo '---- Generate sequence dictionary ----' | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"

java -jar $picard_dir/picard.jar CreateSequenceDictionary O=Genomes/"$VersionMouse"/"$VersionMouse".dict R=Genomes/"$VersionMouse"/"$VersionMouse".fna

echo '---- Generate customized Sanger DB ----' | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"

sh $repository_dir/Preparation_GenerateSangerMouseDB.sh $VersionMouse $temp_dir
rm $temp_dir

echo '---- Generate reference data for msisensor  ----' | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"

msisensor scan -d Genomes/"$VersionMouse"/"$VersionMouse".fna -o Genomes/"$VersionMouse"/"$VersionMouse".microsatellites

echo '---- Generate exons covered by SureSelect ----' | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"
# Download specific .zip-file from https://earray.chem.agilent.com/suredesign/search.htm
# Attention: Most use old versions (mm9/hg19) -> they need liftover
# use "_Regions.bed" for further work - this covers all regions which are targeted in this kit
# Liftover with https://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core using default settings
# Rename to .bed and move to main reference directory
# Included in the data-directory is a version which has alredy been lifted over - nothing more to do but generating the sequence directionary
cp "$repository_dir"/../data/SureSelect_Mouse_All_Exon_V1_mm10.bed Genomes/"$VersionMouse"/SureSelect_Mouse_All_Exon_V1_mm10.bed
java -jar $picard_dir/picard.jar BedToIntervalList I=Genomes/"$VersionMouse"/SureSelect_Mouse_All_Exon_V1_mm10.bed O=Genomes/"$VersionMouse"/SureSelect_Mouse_All_Exon_V1_mm10.bed.list SD=Genomes/"$VersionMouse"/"$VersionMouse".dict

echo '---- Optional for WES: Generating reference data for CopywriteR ----' | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"

Rscript $repository_dir/Preparation_GenerateCopywriterReferences.R

echo '---- Optional for WGS: Generating reference data for HMMCopy ----' | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"

resolution=20000
$hmmcopyutils/util/mappability/generateMap.pl -o Genomes/"$VersionMouse"/"$VersionMouse".fna.map.bw -b -w 150 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y Genomes/"$VersionMouse"/"$VersionMouse".fna
$hmmcopyutils/util/mappability/generateMap.pl -o Genomes/"$VersionMouse"/"$VersionMouse".fna.map.bw -w 150 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y Genomes/"$VersionMouse"/"$VersionMouse".fna
$hmmcopyutils/bin/gcCounter -w $resolution -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y Genomes/"$VersionMouse"/"$VersionMouse".fna > Genomes/"$VersionMouse"/"$VersionMouse".gc.$resolution.wig
$hmmcopyutils/bin/mapCounter -w $resolution -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y Genomes/"$VersionMouse"/"$VersionMouse".fna.map.bw > Genomes/"$VersionMouse"/"$VersionMouse".map.$resolution.wig
rm Genomes/"$VersionMouse"/*.ebwt

echo '---- Finished generating reference data ----' | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "Genomes/"$VersionMouse"/GetReferenceData.txt"

exit 0