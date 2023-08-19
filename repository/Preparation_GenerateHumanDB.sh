#!/bin/bash

##########################################################################################
##
## Preparation_HumanDatabases.sh
##
## Download dbSNP, dbNSFP, COSMIC, gnomAD.
##
##########################################################################################

VersionHuman=$1
config_file=$2

species=Human

#reading configuration from $config_file
. $config_file

wget -nv -P "ref/"$VersionHuman"/" ftp://ftp.ncbi.nih.gov//snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/00-common_all.vcf.gz
wget -nv -P "ref/"$VersionHuman"/" ftp://ftp.ncbi.nih.gov//snp/organisms/human_9606_b150_GRCh38p7/VCF/00-All.vcf.gz
gunzip -f "ref/"$VersionHuman"/00-common_all.vcf.gz"
gunzip -f "ref/"$VersionHuman"/00-All.vcf.gz"
awk '{gsub(/^chr/,""); print}' "ref/"$VersionHuman"/00-common_all.vcf" > "ref/"$VersionHuman"/00-common_all_nochr.vcf"
awk '{gsub(/^chr/,""); print}' "ref/"$VersionHuman"/00-All.vcf" > "ref/"$VersionHuman"/00-All_nochr.vcf"

mv "ref/"$VersionHuman"/00-common_all_nochr.vcf" "ref/"$VersionHuman"/00-common_all.vcf"
bgzip "ref/"$VersionHuman"/00-common_all.vcf"
tabix -p vcf "ref/"$VersionHuman"/00-common_all.vcf.gz"

rm "ref/"$VersionHuman"/00-All.vcf"
mv "ref/"$VersionHuman"/00-All_nochr.vcf" "ref/"$VersionHuman"/00-all.vcf"
bgzip "ref/"$VersionHuman"/00-all.vcf"
tabix -p vcf "ref/"$VersionHuman"/00-all.vcf.gz"

wget -nv -P "ref/"$VersionHuman"/" ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz
wget -nv -P "ref/"$VersionHuman"/" ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz.tbi
wget -nv -P "ref/"$VersionHuman"/" ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz
wget -nv -P "ref/"$VersionHuman"/" ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz.tbi
bcftools norm -m -any "ref/"$VersionHuman"/gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz" -O z -o "ref/"$VersionHuman"/gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.split_multiallelic_sites.vcf.gz"
bcftools norm -m -any "ref/"$VersionHuman"/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz" -O z -o "ref/"$VersionHuman"/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.split_multiallelic_sites.vcf.gz"
mv "ref/"$VersionHuman"/gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.split_multiallelic_sites.vcf.gz" "ref/"$VersionHuman"/gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz"
mv "ref/"$VersionHuman"/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.split_multiallelic_sites.vcf.gz" "ref/"$VersionHuman"/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz"
rm "ref/"$VersionHuman"/gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz.tbi"
rm "ref/"$VersionHuman"/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz.tbi"
tabix -p vcf "ref/"$VersionHuman"/gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz"
tabix -p vcf "ref/"$VersionHuman"/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz"

wget -nv -P "ref/"$VersionHuman"/" https://storage.googleapis.com/gcp-public-data--gnomad/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.vcf.bgz
wget -nv -P "ref/"$VersionHuman"/" https://storage.googleapis.com/gcp-public-data--gnomad/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.vcf.bgz.tbi

wget -nv -P "ref/"$VersionHuman"/" ftp://ftp.ncbi.nlm.nih.gov//pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget -nv -P "ref/"$VersionHuman"/" ftp://ftp.ncbi.nlm.nih.gov//pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
mv "ref/"$VersionHuman"/"clinvar.vcf.gz "ref/"$VersionHuman"/"GRCh38.clinvar.vcf.gz
mv "ref/"$VersionHuman"/"clinvar.vcf.gz.tbi "ref/"$VersionHuman"/"GRCh38.clinvar.vcf.gz.tbi

wget -nv -P "ref/"$VersionHuman"/" ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.basic.annotation.gtf.gz
gunzip -f "ref/"$VersionHuman"/gencode.v28.basic.annotation.gtf.gz"
awk '{gsub(/^chr/,""); print}' "ref/"$VersionHuman"/gencode.v28.basic.annotation.gtf" > "ref/"$VersionHuman"/GRCh38.p12.no_chr.gtf"
cat "ref/"$VersionHuman"/GRCh38.p12.no_chr.gtf" | sort -k1,1 -k4,4n > "ref/"$VersionHuman"/GRCh38.p12.sort.gtf"
mv "ref/"$VersionHuman"/GRCh38.p12.sort.gtf" "ref/"$VersionHuman"/GRCh38.p12.gtf"
bgzip "ref/"$VersionHuman"/GRCh38.p12.gtf"
tabix -p gff "ref/"$VersionHuman"/GRCh38.p12.gtf.gz"

cp $repository_dir"/../data/GRCh38/GRCh38.CosmicCodingMuts.v88.vcf.gz" "ref/"$VersionHuman"/"
cp $repository_dir"/../data/GRCh38/GRCh38.CosmicNonCodingVariants.v88.vcf.gz" "ref/"$VersionHuman"/"
gunzip -f "ref/"$VersionHuman"/GRCh38.CosmicCodingMuts.v88.vcf.gz"
gunzip -f "ref/"$VersionHuman"/GRCh38.CosmicNonCodingVariants.v88.vcf.gz"
sed -i 's/CNT/CNT_NonCoding/g' "ref/"$VersionHuman"/GRCh38.CosmicNonCodingVariants.v88.vcf"
sed -i 's/CNT/CNT_Coding/g' "ref/"$VersionHuman"/GRCh38.CosmicCodingMuts.v88.vcf"
bgzip "ref/"$VersionHuman"/GRCh38.CosmicNonCodingVariants.v88.vcf"
bgzip "ref/"$VersionHuman"/GRCh38.CosmicCodingMuts.v88.vcf"
tabix ref/"$VersionHuman"/GRCh38.CosmicNonCodingVariants.v88.vcf.gz
tabix ref/"$VersionHuman"/GRCh38.CosmicCodingMuts.v88.vcf.gz

echo "Downloading dbNSFP4.4a.zip (~30GB, this can take a while)"
wget -nv -P "ref/"$VersionHuman"/" https://dbnsfp.s3.amazonaws.com/dbNSFP4.4a.zip
cd "ref/"$VersionHuman"/"
unzip dbNSFP4.4a.zip
(zcat dbNSFP4.4a_variant.chr1 | head -n 1; zcat dbNSFP4.4a_variant.chr* | grep -v "^#" ) > dbNSFP4.4a.txt
bgzip dbNSFP4.4a.txt
tabix -s 1 -b 2 -e 2 dbNSFP4.4a.txt.gz
rm dbNSFP4.4a_variant.chr* dbNSFP4.4a.readme.txt dbNSFP4.4_gene.gz dbNSFP4.4_gene.complete.gz
rm LICENSE.txt search_dbNSFP44a.jar search_dbNSFP44a.java search_dbNSFP44a.readme.pdf
rm try.vcf tryhg18.in tryhg19.in dbNSFP4.4a.zip tryhg38.in
cd ../..
