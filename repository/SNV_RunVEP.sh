#!/bin/bash

##########################################################################################
##
## SNV_RunVEP.sh
##
## Running VEP for Mutect2 and Strelka.
##
##########################################################################################

name=$1
config_file=$2
species=$3
input_format=$4
runmode=$5
types=$6

. $config_file

if [ $species = 'Human' ]; then
	species=homo_sapiens
	assembly=GRCh38
	cache_dir=$genome_dir/VEP/
elif [ $species = 'Mouse' ]; then
	species=mus_musculus
	assembly=GRCm38
	cache_dir=$genome_dir/VEP/
fi

if [ $runmode = "MS" ]; then
	types="Tumor Normal"
fi

if [ $input_format = 'Strelka' ]; then
	normalid=NORMAL
	tumorid=TUMOR
elif [ $input_format = 'Mutect2' ]; then
	normalid=Normal
	tumorid=Tumor
elif [ $input_format = 'Manta' ]; then
	normalid=Normal
	tumorid=Tumor
fi

if [ $runmode = 'MS' ]; then

## --polyphen b AND --af_gnomad were removed due to errors!
	$vep_dir/./vep --cache --species $species \
	-i $name/results/$input_format/$name.$input_format.vcf \
	-o $name/results/$input_format/$name.$input_format.vep.vcf \
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


	#perl $vcf2maf_dir/vcf2maf.pl \
	perl $repository_dir/vcf2maf_modified.pl \
	--input-vcf $name/results/$input_format/$name.$input_format.vcf \
	--output-maf $name/results/$input_format/$name.$input_format.vep.maf \
	--vep-path $vep_dir --vep-data $vepdata_dir \
	--ncbi-build $assembly --ref-fasta $genome_file \
	--vcf-normal-id $normalid --vcf-tumor-id $tumorid \
	--tumor-id $name --species $species

	awk 'BEGIN{FS=OFS="\t"} FNR==1{print $0}' $name/results/$input_format/$name.$input_format.vep.maf > $name/results/$input_format/$name.$input_format.vep.tmp
	awk 'BEGIN{FS=OFS="\t"} FNR==2{print $0,$NF="t_maf"}' $name/results/$input_format/$name.$input_format.vep.maf >> $name/results/$input_format/$name.$input_format.vep.tmp
	awk 'BEGIN{FS=OFS="\t"} NR>2{print $0,$NF=$42/$40}' $name/results/$input_format/$name.$input_format.vep.maf >> $name/results/$input_format/$name.$input_format.vep.tmp
	mv $name/results/$input_format/$name.$input_format.vep.tmp $name/results/$input_format/$name.$input_format.vep.maf
fi

if [ $input_format = 'Mutect2' ]; then
	for type in $types;
	do
		(
		## --polyphen b AND --af_gnomad were removed due to errors!
		$vep_dir/./vep --cache --species $species \
		-i $name/results/$input_format/$name.$type.$input_format.vcf \
		-o $name/results/$input_format/$name.$type.$input_format.vep.vcf \
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

		#perl $vcf2maf_dir/vcf2maf.pl \
		perl $repository_dir/vcf2maf_modified.pl \
		--input-vcf $name/results/$input_format/$name.$type.$input_format.vcf \
		--output-maf $name/results/$input_format/$name.$type.$input_format.vep.maf \
		--vep-path $vep_dir --vep-data $vepdata_dir \
		--ncbi-build $assembly --ref-fasta $genome_file \
		--vcf-normal-id Normal --vcf-tumor-id Tumor \
		--tumor-id $name --species $species

		awk 'BEGIN{FS=OFS="\t"} FNR==1{print $0}' $name/results/$input_format/$name.$type.$input_format.vep.maf > $name/results/$input_format/$name.$type.$input_format.vep.maf.tmp
		awk 'BEGIN{FS=OFS="\t"} FNR==2{print $0,$NF="t_maf"}' $name/results/$input_format/$name.$type.$input_format.vep.maf >> $name/results/$input_format/$name.$type.$input_format.vep.maf.tmp

		if [ $type = 'Normal' ]; then
			awk 'BEGIN{FS=OFS="\t"} NR>2{print $0,$NF=0}' $name/results/$input_format/$name.$type.$input_format.vep.maf >> $name/results/$input_format/$name.$type.$input_format.vep.maf.tmp
		fi
		if [ $type = 'Tumor' ]; then
			awk 'BEGIN{FS=OFS="\t"} NR>2{print $0,$NF=$42/$40}' $name/results/$input_format/$name.$type.$input_format.vep.maf >> $name/results/$input_format/$name.$type.$input_format.vep.maf.tmp
		fi
		mv $name/results/$input_format/$name.$type.$input_format.vep.maf.tmp $name/results/$input_format/$name.$type.$input_format.vep.maf
		) &
	done

	wait
fi
