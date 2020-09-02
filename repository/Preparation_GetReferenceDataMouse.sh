#!/bin/bash

##########################################################################################
##
## Preparation_GetReferenceDataMouse.sh
##
## Main routine for the download of all reference data needed for the WES and WGS workflows.
##
##########################################################################################

config_file=$1
temp_dir=$2

species=Mouse

#reading configuration from $config_file
. $config_file

VersionMouse=GRCm38.p6

mkdir -p ref
mkdir -p "ref/"$VersionMouse
mkdir -p "ref/"$VersionMouse/VEP

echo '---- Get reference data ----' | tee "ref/"$VersionMouse"/GetReferenceData.txt"
echo '---- Generate reference data for version '$VersionMouse' ----' | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"

#rerouting STDERR to report file
exec 2>> "ref/"$VersionMouse"/GetReferenceData.txt"

echo '---- Copying over files from repository ----' | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"

cp $repository_dir"/../data/GRCm38/GRCm38.canonical_chromosomes.bed" "ref/"$VersionMouse"/"
bgzip "ref/"$VersionMouse"/GRCm38/GRCm38.canonical_chromosomes.bed"
tabix -p bed "ref/"$VersionMouse"/GRCm38/GRCm38.canonical_chromosomes.bed.gz"

cp $repository_dir"/../data/GRCm38/GRCm38.Census_allMon_Jan_15_11_46_18_2018_mouse.tsv" "ref/"$VersionMouse"/"

cp $repository_dir"/../data/GRCm38/GRCm38.bammatcher_docker.conf" "ref/"$VersionMouse"/"
cp $repository_dir"/../data/GRCm38/GRCm38.bammatcher_bash.conf" "ref/"$VersionMouse"/"

cp $repository_dir"/../data/GRCm38/GRCm38.AgilentProbeGaps.txt" "ref/"$VersionMouse"/"
cp $repository_dir"/../data/GRCm38/GRCm38.Genecode_M20_Exons.rds" "ref/"$VersionMouse"/"
cp $repository_dir"/../data/GRCm38/GRCm38.Genecode_M20_Genes.rds" "ref/"$VersionMouse"/"

cp $repository_dir"/../data/GRCm38/GRCm38.RefFlat" "ref/"$VersionMouse"/"

cp $repository_dir"/../data/Samples.tsv" "ref/"$VersionMouse"/"

cp $repository_dir"/../data/GRCm38/GRCm38.SureSelect_Mouse_All_Exon_V1.bed" ref/"$VersionMouse"/

echo '---- Downloading reference genome ----' | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"

wget -nv -P "ref/"$VersionMouse "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.8_"$VersionMouse"/GCA_000001635.8_"$VersionMouse"_genomic.fna.gz"
gunzip "ref/"$VersionMouse"/GCA_000001635.8_"$VersionMouse"_genomic.fna.gz"
mv "ref/"$VersionMouse"/GCA_000001635.8_"$VersionMouse"_genomic.fna" "ref/"$VersionMouse"/"$VersionMouse".fna"

echo '---- Generate BWA Index ----' | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"

sh $repository_dir/Preparation_GenerateBWAIndex.sh $VersionMouse $config_file $species

echo '---- Generate sequence dictionary ----' | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"

java -Dpicard.useLegacyParser=false -jar $picard_dir/picard.jar CreateSequenceDictionary \
-O "ref/"$VersionMouse"/"$VersionMouse".dict" \
-R "ref/"$VersionMouse"/"$VersionMouse".fna"

echo '---- Generate exons covered by SureSelect ----' | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"
# Download specific .zip-file from https://earray.chem.agilent.com/suredesign/search.htm
# Attention: Most use old versions (mm9/hg19) -> they need liftover
# use "_Regions.bed" for further work - this covers all regions which are targeted in this kit
# Liftover with https://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core using default settings
# Rename to .bed and move to main reference directory
# Included in the data-directory is a version which has already been lifted over - nothing more to do but generating the sequence directionary
java -Dpicard.useLegacyParser=false -jar $picard_dir/picard.jar BedToIntervalList \
-I "ref/"$VersionMouse"/GRCm38.SureSelect_Mouse_All_Exon_V1.bed" \
-O "ref/"$VersionMouse"/GRCm38.SureSelect_Mouse_All_Exon_V1.bed.list" \
-SD "ref/"$VersionMouse"/"$VersionMouse".dict"

echo '---- Downloading reference genome (for VEP) ----' | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"

wget -nv -P "ref/"$VersionMouse"/VEP" ftp://ftp.ensembl.org/pub/release-96/variation/indexed_vep_cache/mus_musculus_vep_96_GRCm38.tar.gz
tar -xzf "ref/"$VersionMouse"/VEP/mus_musculus_vep_96_GRCm38.tar.gz"
mv mus_musculus "ref/"$VersionMouse"/VEP/"
rm "ref/"$VersionMouse"/VEP/mus_musculus_vep_96_GRCm38.tar.gz"

echo '---- Generate customized Sanger DB ----' | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"

sh $repository_dir/Preparation_GenerateSangerMouseDB.sh $VersionMouse $temp_dir
rm -rf $temp_dir

wget -nv -c -r -P ref/"$VersionMouse"/ ftp://ftp-mouse.sanger.ac.uk/REL-1807-SNPs_Indels/mgp.v6.merged.norm.snp.indels.sfiltered.vcf.gz
wget -nv -c -r -P ref/"$VersionMouse"/ ftp://ftp-mouse.sanger.ac.uk/REL-1807-SNPs_Indels/mgp.v6.merged.norm.snp.indels.sfiltered.vcf.gz.tbi
mv ref/"$VersionMouse"/ftp-mouse.sanger.ac.uk/REL-1807-SNPs_Indels/mgp.v6.merged.norm.snp.indels.sfiltered.vcf.gz ref/"$VersionMouse"/MGP.v6.snp_and_indels.vcf.gz
mv ref/"$VersionMouse"/ftp-mouse.sanger.ac.uk/REL-1807-SNPs_Indels/mgp.v6.merged.norm.snp.indels.sfiltered.vcf.gz.tbi ref/"$VersionMouse"/MGP.v6.snp_and_indels.vcf.gz.tbi
rm -r ref/"$VersionMouse"/ftp-mouse.sanger.ac.uk/REL-1807-SNPs_Indels/
rm -r ref/"$VersionMouse"/ftp-mouse.sanger.ac.uk/

bcftools view -s ^CAST_EiJ,SPRET_EiJ,PWK_PhJ,WSB_EiJ,MOLF_EiJ,ZALENDE_EiJ,LEWES_EiJ --min-ac=1 --no-update ref/"$VersionMouse"/MGP.v6.snp_and_indels.vcf.gz -O z -o ref/"$VersionMouse"/MGP.v6.snp_and_indels.exclude_wild.vcf.gz
tabix -p vcf ref/"$VersionMouse"/MGP.v6.snp_and_indels.exclude_wild.vcf.gz

rm "ref/"$VersionMouse"/MGP.v6.snp_and_indels.vcf.gz"
rm "ref/"$VersionMouse"/MGP.v6.snp_and_indels.vcf.gz.tbi"

echo '---- Generate reference data for msisensor  ----' | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"

msisensor scan -d "ref/"$VersionMouse"/"$VersionMouse".fna" -o "ref/"$VersionMouse"/"$VersionMouse".microsatellites"

echo '---- Optional for WES: Generating reference data for CopywriteR ----' | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"

Rscript $repository_dir/Preparation_GenerateCopywriterReferences.R $VersionMouse

echo '---- Optional for WGS: Generating reference data for HMMCopy ----' | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"

$hmmcopyutils_dir/util/mappability/generateMap.pl -o "ref/"$VersionMouse"/"$VersionMouse".fna.map.bw" -b -w 150 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y "ref/"$VersionMouse"/"$VersionMouse".fna"
$hmmcopyutils_dir/util/mappability/generateMap.pl -o "ref/"$VersionMouse"/"$VersionMouse".fna.map.bw" -w 150 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y "ref/"$VersionMouse"/"$VersionMouse".fna"

for resolution in 1000 10000 20000 50000;
do
	$hmmcopyutils_dir/bin/gcCounter -w $resolution -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y "ref/"$VersionMouse"/"$VersionMouse".fna" > "ref/"$VersionMouse"/"$VersionMouse".gc.$resolution.wig"
	$hmmcopyutils_dir/bin/mapCounter -w $resolution -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y "ref/"$VersionMouse"/"$VersionMouse".fna.map.bw" > "ref/"$VersionMouse"/"$VersionMouse".map.$resolution.wig"
done
rm "ref/"$VersionMouse"/*.ebwt"

echo '---- Finished generating reference data ----' | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"
date | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"
echo 'DONE' | tee -a "ref/"$VersionMouse"/GetReferenceData.txt"

exit 0
