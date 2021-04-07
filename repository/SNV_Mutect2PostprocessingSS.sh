#!/bin/bash

##########################################################################################
##
## SNV_Mutect2Postprocessing.sh
##
## Postprocessing for Mutect2 in single-sample mode.
##
##########################################################################################

name=$1
species=$2
config_file=$3
filtering=$4
artefact_type=$5
GATK=$6
type=$7

. $config_file

echo '---- Mutect2 SS Postprocessing I (OrientationFilter, Indel size selection, filtering) ----' | tee -a ${name}/results/QC/${name}.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a ${name}/results/QC/${name}.report.txt


# this will change the "PASS" flag, which is filtered later with SnpSift.jar filter
if [ $artefact_type = 'no' ]; then
	# just copy without filtering
	cp ${name}/results/Mutect2/${name}.${type}.m2.vcf ${name}/results/Mutect2/${name}.${type}.m2.filt.vcf

elif [ $artefact_type = 'yes' ]; then
	# We expect the ob-file to already be available
	# filter with ob-priors
	java -jar $GATK_dir/gatk.jar FilterMutectCalls \
	--variant ${name}/results/Mutect2/${name}.${type}.m2.vcf \
	--output ${name}/results/Mutect2/${name}.${type}.m2.filt.vcf \
	--reference $genome_file \
	-ob-priors ${name}/results/Mutect2/${name}.${type}.m2.read-orientation-model.tar.gz
fi

# output filtering statistics
grep "^[^#;]" $name/results/Mutect2/$name.$type.m2.filt.vcf | cut -f 7 | sort | uniq -c | sort -nr > $name/results/Mutect2/$name.$type.m2.filt.filtersummary.txt

java -jar $snpeff_dir/SnpSift.jar extractFields \
$name/results/Mutect2/"$name".$type.m2.filt.vcf \
CHROM POS REF ALT "GEN["$type"].AF" "GEN["$type"].AD[0]" \
"GEN["$type"].AD[1]" MMQ[1] MBQ[1] \
> $name/results/Mutect2/$name.$type.Mutect2.Positions.txt

java -jar $GATK_dir/gatk.jar SelectVariants --max-indel-size 10 \
-V $name/results/Mutect2/$name.$type.m2.filt.vcf \
-output $name/results/Mutect2/$name.$type.m2.filt.selected.vcf

if [ $filtering = 'soft' ]; then
	cat $name/results/Mutect2/$name.$type.m2.filt.selected.vcf \
	| java -jar $snpeff_dir/SnpSift.jar filter \
	"( ( FILTER = 'PASS') & (GEN[$type].AF >= 0.05) & \
	(GEN[$type].AD[1] >= 2) & (GEN[$type].AD[0] + GEN[$type].AD[1] >= 5) )" \
	> $name/results/Mutect2/$name.$type.m2.postprocessed.vcf
elif [ $filtering = 'hard' ]; then
	cat $name/results/Mutect2/$name.$type.m2.filt.selected.vcf \
	| java -jar $snpeff_dir/SnpSift.jar filter \
	"( ( FILTER = 'PASS') & (GEN[$type].AF >= 0.1) & \
	(GEN[$type].AD[1] >= 2) & (GEN[$type].AD[0] + GEN[$type].AD[1] >= 10) )" \
	> $name/results/Mutect2/$name.$type.m2.postprocessed.vcf
elif [ $filtering = 'none' ]; then
	cat $name/results/Mutect2/$name.$type.m2.filt.selected.vcf \
	| java -jar $snpeff_dir/SnpSift.jar filter \
	"( ( FILTER = 'PASS' ) )" \
	> $name/results/Mutect2/$name.$type.m2.postprocessed.vcf
fi

echo '---- Mutect2 SS Postprocessing II (Filtering out known SNV/Indel using dbSNP or the Sanger Mouse database) ----' | tee -a $name/results/QC/$name.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

bgzip -f $name/results/Mutect2/$name.$type.m2.postprocessed.vcf

tabix -p vcf $name/results/Mutect2/$name.$type.m2.postprocessed.vcf.gz

if [ $filtering = 'soft' ]; then
	bcftools isec -C -c none -O z -w 1 \
	-o $name/results/Mutect2/$name.$type.m2.postprocessed.snp_removed.vcf.gz \
	$name/results/Mutect2/$name.$type.m2.postprocessed.vcf.gz \
	$snp_file
elif [ $filtering = 'hard' ]; then
	bcftools isec -C -c none -O z -w 1 \
	-o $name/results/Mutect2/$name.$type.m2.postprocessed.snp_removed.vcf.gz \
	$name/results/Mutect2/$name.$type.m2.postprocessed.vcf.gz \
	$alternate_snp_file
elif [ $filtering = 'none' ]; then
	cp $name/results/Mutect2/$name.$type.m2.postprocessed.vcf.gz \
	$name/results/Mutect2/$name.$type.m2.postprocessed.snp_removed.vcf.gz
fi

bcftools norm -m -any $name/results/Mutect2/$name.$type.m2.postprocessed.snp_removed.vcf.gz -O z -o $name/results/Mutect2/$name.$type.Mutect2.vcf.gz

gunzip -f $name/results/Mutect2/$name.$type.Mutect2.vcf.gz

echo '---- Mutect2 SS Postprocessing III (Annotate calls) ----' | tee -a $name/results/QC/$name.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

if [ $species = 'Human' ]; then

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar annotate \
	$dbsnp_file $name/results/Mutect2/$name.$type.Mutect2.vcf \
	> $name/results/Mutect2/$name.$type.Mutect2.ann1.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar \
	annotate $cosmiccoding_file $name/results/Mutect2/$name.$type.Mutect2.ann1.vcf \
	> $name/results/Mutect2/$name.$type.Mutect2.ann2.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar  \
	annotate $cosmicnoncoding_file $name/results/Mutect2/$name.$type.Mutect2.ann2.vcf \
	> $name/results/Mutect2/$name.$type.Mutect2.ann3.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar  \
	annotate $clinvar_file $name/results/Mutect2/$name.$type.Mutect2.ann3.vcf \
	> $name/results/Mutect2/$name.$type.Mutect2.ann4.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar  \
	annotate $gnomadexome_file $name/results/Mutect2/$name.$type.Mutect2.ann4.vcf \
	> $name/results/Mutect2/$name.$type.Mutect2.ann5.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar  \
	annotate $gnomadgenome_file $name/results/Mutect2/$name.$type.Mutect2.ann5.vcf \
	> $name/results/Mutect2/$name.$type.Mutect2.ann6.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar \
	DbNSFP -db $dbnsfp_file -v $name/results/Mutect2/"$name".$type.Mutect2.ann6.vcf \
	-f MetaLR_pred,MetaSVM_pred,SIFT_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,PROVEAN_pred \
	> $name/results/Mutect2/"$name".$type.Mutect2.ann7.vcf

	java -Xmx16g -jar $snpeff_dir/snpEff.jar $snpeff_version -canon \
	-csvStats $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf.stats \
	$name/results/Mutect2/"$name".$type.Mutect2.ann7.vcf \
	> $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf

	cat $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf \
	| $snpeff_dir/scripts/vcfEffOnePerLine.pl \
	> $name/results/Mutect2/$name.$type.Mutect2.annotated.one.vcf

	java -jar $snpeff_dir/SnpSift.jar extractFields \
	$name/results/Mutect2/$name.$type.Mutect2.annotated.one.vcf \
	 CHROM POS REF ALT "GEN[$type].AF" "GEN[$type].AD[0]" \
	 "GEN[$type].AD[1]" ANN[*].GENE  ANN[*].EFFECT ANN[*].IMPACT \
	 ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
	 dbNSFP_MetaLR_pred dbNSFP_MetaSVM_pred ID CAF G5 AC AN AF CNT_Coding \
	 CNT_NonCoding CLNDN CLNSIG CLNREVSTAT dbNSFP_SIFT_pred \
	 dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_pred dbNSFP_PROVEAN_pred \
	 > $name/results/Mutect2/$name.$type.Mutect2.txt

elif [ $species = 'Mouse' ]; then

	java -Xmx16g -jar $snpeff_dir/snpEff.jar $snpeff_version -canon \
	-csvStats $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf.stats \
	$name/results/Mutect2/$name.$type.Mutect2.vcf \
	> $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf

	cat $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf \
	| $snpeff_dir/scripts/vcfEffOnePerLine.pl \
	> $name/results/Mutect2/$name.$type.Mutect2.annotated.one.vcf

	java -jar $snpeff_dir/SnpSift.jar extractFields \
	$name/results/Mutect2/$name.$type.Mutect2.annotated.one.vcf \
	CHROM POS REF ALT "GEN[$type].AF" "GEN[$type].AD[0]" \
	"GEN[$type].AD[1]" ANN[*].GENE  ANN[*].EFFECT ANN[*].IMPACT \
	ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
	> $name/results/Mutect2/$name.$type.Mutect2.txt

fi

sh $repository_dir/SNV_CleanUp.sh $name SS $type $species
