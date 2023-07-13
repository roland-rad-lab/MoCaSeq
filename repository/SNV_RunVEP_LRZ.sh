#!/bin/bash

##########################################################################################
##
## SNV_RunVEP.sh
##
## Running VEP for Mutect2 and Strelka.
##
##########################################################################################
name=PCSI_0410_Ag_M_526
config_file=/opt/MoCaSeq/config.sh
species=Human
input_format='Mutect2'
runmode='MS'
types="Tumor Normal"

name=$1
config_file=$2
species=$3
input_format=$4
runmode=$5
types=$6

. $config_file

if [ $species = 'Human' ]; then
	species_tax=homo_sapiens
	assembly=GRCh38
	cache_dir=$genome_dir/VEP/
elif [ $species = 'Mouse' ]; then
	species_tax=mus_musculus
	assembly=GRCm38
	cache_dir=$genome_dir/VEP/
fi

if [ $runmode = "MS" ]; then
	types="Tumor Normal"
fi

if [ ${input_format} = 'Strelka' ]; then
	normalid=NORMAL
	tumorid=TUMOR
elif [ ${input_format} = 'Mutect2' ]; then
	normalid=Normal
	tumorid=Tumor
elif [ ${input_format} = 'Manta' ]; then
	normalid=Normal
	tumorid=Tumor
fi

if [ $runmode = 'MS' ]; then

	# this script will extract the transcripts identified as canonical from the main (unfiltered) mutation file
	# this will be used to force vcf2maf to convert those instead of internally defined ones,
	# so we are consistent between TXT and MAF
	# output file: NAME/results/Mutect2/NAME.Mutect2.canonicalTranscripts.txt
	Rscript $repository_dir/SNV_VCF-TXT_to_CanonicalTranscripts.R \
	${name}/results/${input_format}/${name}.${input_format}.txt \
	$species

	# add a custom header to the vcf file to facilitate the merging of VCF and MAF
	python $repository_dir/SNV_Add_ID_To_VCF.py --i "${name}/results/${input_format}/${name}.${input_format}.vcf" --o "${name}/results/${input_format}/${name}.${input_format}.header.vcf"

	# --polyphen b AND --af_gnomad were removed due to errors!
	$vep_dir/./vep --cache --species $species_tax \
	-i ${name}/results/${input_format}/${name}.${input_format}.header.vcf \
	-o ${name}/results/${input_format}/${name}.${input_format}.vep.vcf \
	--fasta $genome_file --assembly $assembly \
	--offline --no_progress --no_stats \
	--buffer_size 5000 --sift b --ccds --uniprot --hgvs \
	--symbol --numbers --domains --gene_phenotype --canonical \
	--protein --biotype --uniprot --tsl --pubmed --variant_class \
	--shift_hgvs 1 --check_existing --total_length --allele_number \
	--no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele \
	--pick_order canonical,tsl,biotype,rank,ccds,length --format vcf \
	--fork 4 --cache_version 96 --af --af_1kg --af_esp \
	--force_overwrite --dir_cache $cache_dir

	# this is VEP (as above) + MAF, therefore the input is the classic VCF
	perl $repository_dir/vcf2maf_modified.pl \
	--input-vcf ${name}/results/${input_format}/${name}.${input_format}.vcf \
	--output-maf ${name}/results/${input_format}/${name}.${input_format}.vep.maf \
	--vep-path $vep_dir --vep-data $vepdata_dir \
	--ncbi-build $assembly --ref-fasta $genome_file \
	--vcf-normal-id $normalid --vcf-tumor-id $tumorid \
	--tumor-id ${name} --species $species_tax \
	--custom-enst ${name}/results/${input_format}/${name}.${input_format}.canonicalTranscripts.txt \
	 --retain-fmt MERGEID --retain-info MERGEID

	 # remove the (buggy and undefined) output
	 rm -f ${name}/results/${input_format}/${name}.${input_format}.vep.vep.maf
	 rm -f ${name}/results/${input_format}/${name}.${input_format}.vep.vep.vcf

	# 1. extract the version commend line
	awk 'BEGIN{FS=OFS="\t"} FNR==1{print $0}' ${name}/results/${input_format}/${name}.${input_format}.vep.maf > ${name}/results/${input_format}/${name}.${input_format}.vep.tmp
	# 2. add t_maf as column
	awk 'BEGIN{FS=OFS="\t"} FNR==2{print $0,$NF="t_maf"}' ${name}/results/${input_format}/${name}.${input_format}.vep.maf >> ${name}/results/${input_format}/${name}.${input_format}.vep.tmp
	# 3. calculate the frequency based on column 42/40
	awk 'BEGIN{FS=OFS="\t"} NR>2{print $0,$NF=$42/$40}' ${name}/results/${input_format}/${name}.${input_format}.vep.maf >> ${name}/results/${input_format}/${name}.${input_format}.vep.tmp

	# 4. merge everything back together
	mv ${name}/results/${input_format}/${name}.${input_format}.vep.tmp ${name}/results/${input_format}/${name}.${input_format}.vep.maf
fi

if [ ${input_format} = 'Mutect2' ]; then
	for type in $types;
	do
		(

		# this script will extract the transcripts identified as canonical from the main (unfiltered) mutation file
		# this will be used to force vcf2maf to convert those instead of internally defined ones,
		# so we are consistent between TXT and MAF
		# output file: NAME/results/Mutect2/NAME.TYPE.Mutect2.canonicalTranscripts.txt
		Rscript $repository_dir/SNV_VCF-TXT_to_CanonicalTranscripts.R ${name}/results/${input_format}/${name}.$type.${input_format}.txt $species

		# add a custom header to the vcf file to facilitate the merging of VCF and MAF
		python $repository_dir/SNV_Add_ID_To_VCF.py --i "${name}/results/${input_format}/${name}.$type.${input_format}.vcf" --o "${name}/results/${input_format}/${name}.$type.${input_format}.mergeid.vcf"

		## --polyphen b AND --af_gnomad were removed due to errors!
		$vep_dir/./vep --cache --species $species_tax \
		-i ${name}/results/${input_format}/${name}.$type.${input_format}.mergeid.vcf \
		-o ${name}/results/${input_format}/${name}.$type.${input_format}.vep.vcf \
		--fasta $genome_file --assembly $assembly \
		--offline --no_progress --no_stats \
		--buffer_size 5000 --sift b --ccds --uniprot --hgvs \
		--symbol --numbers --domains --gene_phenotype --canonical \
		--protein --biotype --uniprot --tsl --pubmed --variant_class \
		--shift_hgvs 1 --check_existing --total_length --allele_number \
		--no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele \
		--pick_order canonical,tsl,biotype,rank,ccds,length --format vcf \
		--fork 4 --cache_version 96 --af --af_1kg --af_esp \
		--force_overwrite --dir_cache $cache_dir

		# this is VEP (as above) + MAF, therefore the input is the classic VCF
		perl $repository_dir/vcf2maf_modified.pl \
		--input-vcf ${name}/results/${input_format}/${name}.$type.${input_format}.vcf \
		--output-maf ${name}/results/${input_format}/${name}.$type.${input_format}.vep.maf \
		--vep-path $vep_dir --vep-data $vepdata_dir \
		--ncbi-build $assembly --ref-fasta $genome_file \
		--vcf-normal-id Normal --vcf-tumor-id Tumor \
		--tumor-id ${name} --species $species_tax \
		--custom-enst ${name}/results/${input_format}/${name}.$type.${input_format}.canonicalTranscripts.txt \
		 --retain-fmt MERGEID --retain-info MERGEID

		 # remove the (buggy and undefined) output
		 rm -f ${name}/results/${input_format}/${name}.${input_format}.vep.vep.maf
		 rm -f ${name}/results/${input_format}/${name}.${input_format}.vep.vep.vcf

		 # 1. extract the version commend line
		awk 'BEGIN{FS=OFS="\t"} FNR==1{print $0}' ${name}/results/${input_format}/${name}.$type.${input_format}.vep.maf > ${name}/results/${input_format}/${name}.$type.${input_format}.vep.maf.tmp
		# 2. add t_maf as column
		awk 'BEGIN{FS=OFS="\t"} FNR==2{print $0,$NF="t_maf"}' ${name}/results/${input_format}/${name}.$type.${input_format}.vep.maf >> ${name}/results/${input_format}/${name}.$type.${input_format}.vep.maf.tmp

		 # 3. calculate the frequency based on column 42/40
		if [ $type = 'Normal' ]; then
			awk 'BEGIN{FS=OFS="\t"} NR>2{print $0,$NF=0}' ${name}/results/${input_format}/${name}.$type.${input_format}.vep.maf >> ${name}/results/${input_format}/${name}.$type.${input_format}.vep.maf.tmp
		fi
		if [ $type = 'Tumor' ]; then
			awk 'BEGIN{FS=OFS="\t"} NR>2{print $0,$NF=$42/$40}' ${name}/results/${input_format}/${name}.$type.${input_format}.vep.maf >> ${name}/results/${input_format}/${name}.$type.${input_format}.vep.maf.tmp
		fi
		# 4. merge everything back together
		mv ${name}/results/${input_format}/${name}.$type.${input_format}.vep.maf.tmp ${name}/results/${input_format}/${name}.$type.${input_format}.vep.maf
		) &
	done

	wait
fi
