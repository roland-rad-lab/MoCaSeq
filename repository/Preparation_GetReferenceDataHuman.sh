#!/bin/bash

##########################################################################################
##
## Preparation_GetReferenceDataHuman.sh
##
## Main routine for the download of all reference data needed for the WES and WGS workflows.
##
##########################################################################################

#config_file=$1
#temp_dir=$2
config_file=/opt/MoCaSeq/config.sh
temp_dir=/var/pipeline/temp

species=Human

#reading configuration from $config_file
. $config_file

VersionHuman=GRCh38.p12

mkdir -p ref
mkdir -p "ref/"$VersionHuman
mkdir -p "ref/"$VersionHuman/VEP

echo '---- Get reference data ----' | tee "ref/"$VersionHuman"/GetReferenceData.txt"
echo '---- Generate reference data for version '$VersionHuman' ----' | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"
date | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"

#rerouting STDERR to report file
exec 2>> "ref/"$VersionHuman"/GetReferenceData.txt"

echo '---- Copying over files from repository ----' | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"
date | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"

cp $repository_dir"/../data/GRCh38/GRCh38.canonical_chromosomes.bed" "ref/"$VersionHuman"/"
bgzip "ref/"$VersionHuman"/GRCh38.canonical_chromosomes.bed"
tabix -p bed "ref/"$VersionHuman"/GRCh38.canonical_chromosomes.bed.gz"

cp $repository_dir"/../data/GRCh38/GRCh38.Census_allMon_Feb_11_14_43_15_2019.tsv" "ref/"$VersionHuman"/"

cp $repository_dir"/../data/GRCh38/GRCh38.bammatcher_docker.conf" "ref/"$VersionHuman"/"
cp $repository_dir"/../data/GRCh38/GRCh38.bammatcher_bash.conf" "ref/"$VersionHuman"/"

cp $repository_dir"/../data/GRCh38/gencode_humanV31_exons.rds" "ref/"$VersionHuman"/"
cp $repository_dir"/../data/GRCh38/gencode_humanV31_genes.rds" "ref/"$VersionHuman"/"

cp $repository_dir"/../data/GRCh38/GRCh38.RefFlat" "ref/"$VersionHuman"/"

cp $repository_dir"/../data/Samples.tsv" "ref/"$VersionHuman"/"

cp $repository_dir"/../data/GRCh38/GRCh38.TruSight_Panel_V2_Gene_List.txt" "ref/"$VersionHuman"/"

cp $repository_dir"/../data/GRCh38/GRCh38.SureSelect_Human_All_Exon_V5_hg38.bed" "ref/$VersionHuman"/
cp $repository_dir"/../data/GRCh38/GRCh38_SureSelect.S31285117_MergedProbes_All_Exons_V7.bed" "ref/$VersionHuman"/


echo '---- Downloading centromere files ----' | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"
date | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"

wget -nv -P "ref/"$VersionHuman"/" "https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt"
awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3}' "ref/"$VersionHuman"/GRCh38.GCA_000001405.2_centromere_acen.txt" > "ref/"$VersionHuman"/GRCh38.centromers.txt"
rm "ref/"$VersionHuman"/GRCh38.GCA_000001405.2_centromere_acen.txt"
sed -i 's/chr//g' "ref/"$VersionHuman"/GRCh38.centromers.txt"
sed -i 's/Chr/chr/g' "ref/"$VersionHuman"/GRCh38.centromers.txt"
sed -i 's/Start/start/g' "ref/"$VersionHuman"/GRCh38.centromers.txt"
sed -i 's/End/stop/g' "ref/"$VersionHuman"/GRCh38.centromers.txt"




echo '---- Downloading reference genome ----' | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"
date | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"

# FINISHED BUT NOT TESTED
cd "ref/"$VersionHuman

/opt/bwa-0.7.17/bwakit/run-gen-ref hs38DH

wget -nv ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz

gunzip -f GCF_000819615.1_ViralProj14015_genomic.fna.gz

sed '1c>Phix' GCF_000819615.1_ViralProj14015_genomic.fna > GCF_000819615.1_ViralProj14015_genomic.temp.fna

cat hs38DH.fa GCF_000819615.1_ViralProj14015_genomic.temp.fna > temp.fa

cat temp.fa | sed 's/>chr/>/g' > $VersionHuman'.fna'

sed -i 's/chr//g' hs38DH.fa.alt

mv hs38DH.fa.alt $VersionHuman'.fna.alt'

samtools faidx $VersionHuman'.fna'

rm GCF_000819615.1_ViralProj14015_genomic.fna GCF_000819615.1_ViralProj14015_genomic.temp.fna hs38DH.fa temp.fa

cd ../..






echo '---- Generate BWA Index (~2-4h) ----' | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"
date | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"

# FINISHED AND TESTED
sh $repository_dir/Preparation_GenerateBWAIndex.sh $VersionHuman $config_file $species



echo '---- Generate sequence dictionary ----' | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"
date | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"

FINISHED AND TESTED
java -Dpicard.useLegacyParser=false -jar $picard_dir/picard.jar CreateSequenceDictionary \
-O "ref/"$VersionHuman"/"$VersionHuman".dict" \
-R "ref/"$VersionHuman"/"$VersionHuman".fna"



echo '---- Generate exons covered by SureSelect ----' | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"
date | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"

java -Dpicard.useLegacyParser=false -jar $picard_dir/picard.jar BedToIntervalList \
-I "ref/"$VersionHuman"/GRCh38.SureSelect_Human_All_Exon_V5_hg38.bed" \
-O "ref/"$VersionHuman"/GRCh38.SureSelect_Human_All_Exon_V5_hg38.bed.list" \
-SD "ref/"$VersionHuman"/"$VersionHuman".dict"

java -Dpicard.useLegacyParser=false -jar $picard_dir/picard.jar BedToIntervalList \
-I "ref/"$VersionHuman"/GRCh38_SureSelect.S31285117_MergedProbes_All_Exons_V7.bed" \
-O "ref/"$VersionHuman"/GRCh38_SureSelect.S31285117_MergedProbes_All_Exons_V7.bed.list" \
-SD "ref/"$VersionHuman"/"$VersionHuman".dict"


echo '---- Downloading reference genome (for VEP) ----' | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"
date | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"

# FINISHED AND TESTED
wget -nv -P "ref/"$VersionHuman"/VEP" ftp://ftp.ensembl.org/pub/release-96/variation/indexed_vep_cache/homo_sapiens_vep_96_GRCh38.tar.gz
tar -xzf "ref/"$VersionHuman"/VEP/homo_sapiens_vep_96_GRCh38.tar.gz"
mv homo_sapiens "ref/"$VersionHuman"/VEP/"
rm "ref/"$VersionHuman"/VEP/mus_musculus_vep_96_GRCm38.tar.gz"



echo '---- Generating multiple database files (dbSNP, dbNSFP, COSMIC, gnomAD) ----' | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"
date | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"

sh $repository_dir/Preparation_GenerateHumanDB.sh $VersionHuman $config_file

echo '---- Finished generating reference data ----' | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"
date | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"
echo 'DONE' | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"





echo '---- Generate reference data for msisensor  ----' | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"
date | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"

msisensor scan -d "ref/"$VersionHuman"/"$VersionHuman".fna" -o "ref/"$VersionHuman"/"$VersionHuman".microsatellites"



echo '---- Optional for WES: Generating reference data for CopywriteR ----' | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"
date | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"

Rscript $repository_dir/Preparation_GenerateCopywriterReferences.R $VersionHuman



echo '---- Optional for WGS: Generating reference data for HMMCopy ----' | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"
date | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"

$hmmcopyutils_dir/util/mappability/generateMap.pl -o "ref/"$VersionHuman"/"$VersionHuman".fna.map.bw" -b -w 150 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y "ref/"$VersionHuman"/"$VersionHuman".fna"
$hmmcopyutils_dir/util/mappability/generateMap.pl -o "ref/"$VersionHuman"/"$VersionHuman".fna.map.bw" -w 150 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y "ref/"$VersionHuman"/"$VersionHuman".fna"

for resolution in 1000 10000 20000 50000;
do
	$hmmcopyutils_dir/bin/gcCounter -w $resolution -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y "ref/"$VersionHuman"/"$VersionHuman".fna" > "ref/"$VersionHuman"/"$VersionHuman".gc.$resolution.wig"
	$hmmcopyutils_dir/bin/mapCounter -w $resolution -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y "ref/"$VersionHuman"/"$VersionHuman".fna.map.bw" > "ref/"$VersionHuman"/"$VersionHuman".map.$resolution.wig"
done
rm "ref/"$VersionHuman"/*.ebwt"

echo '---- Finished generating reference data ----' | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"
date | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"
echo 'DONE' | tee -a "ref/"$VersionHuman"/GetReferenceData.txt"

exit 0
