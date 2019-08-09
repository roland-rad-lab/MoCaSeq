#!/bin/bash

##########################################################################################
##
## SNV_Mutect2Postprocessing.sh
##
## Postprocessing for Mutect2.
##
##########################################################################################

name=$1
species=$2
config_file=$3
filtering=$4
artefact_type=$5
GATK=$6

. $config_file

echo '---- Mutect2 Postprocessing I (OrientationFilter, Indel size selection, filtering) ----' | tee -a $name/results/QC/$name.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

java -jar $GATK_dir/gatk.jar FilterMutectCalls \
--variant $name/results/Mutect2/"$name".m2.vcf \
--output $name/results/Mutect2/"$name".m2.filt.vcf \
--reference $genome_file

if [ $artefact_type = 'none' ]; then
	cp $name/results/Mutect2/$name.m2.filt.vcf \
	$name/results/Mutect2/$name.m2.filt.AM.vcf
elif [ $artefact_type = 'CT' ]; then
	java -jar $GATK_dir/gatk.jar FilterByOrientationBias \
	-V $name/results/Mutect2/$name.m2.filt.vcf -P \
	$name/results/QC/$name.Tumor.bam.artifacts.pre_adapter_detail_metrics \
	--artifact-modes C/T --output $name/results/Mutect2/$name.m2.filt.AM.vcf
elif [ $artefact_type = 'GT' ]; then
	java -jar $GATK_dir/gatk.jar FilterByOrientationBias \
	-V $name/results/Mutect2/$name.m2.filt.vcf -P \
	$name/results/QC/$name.Tumor.bam.artifacts.pre_adapter_detail_metrics \
	--artifact-modes G/T --output $name/results/Mutect2/$name.m2.filt.AM.vcf
fi

if [ $artefact_type = 'none' ]; then
	cp $name/results/Mutect2/$name.m2.filt.AM.vcf \
	$name/results/Mutect2/$name.m2.filt.AM.filtered.vcf
elif [ $artefact_type = 'CT' ] || [ $artefact_type = 'GT' ]; then
	cat $name/results/Mutect2/$name.m2.filt.AM.vcf \
	| java -jar $snpeff_dir/SnpSift.jar filter \
	"( ( ( FILTER = 'PASS'  ) & (exists GEN[Tumor].OBP) & \
	(GEN[Tumor].OBP <= 0.05) ) | ( ( FILTER = 'PASS'  ) ) )" \
	> $name/results/Mutect2/$name.m2.filt.AM.filtered.vcf
fi

java -jar $GATK_dir/gatk.jar SelectVariants --max-indel-size 10 \
-V $name/results/Mutect2/$name.m2.filt.AM.filtered.vcf \
-output $name/results/Mutect2/$name.m2.filt.AM.filtered.selected.vcf

if [ $filtering = 'all' ]; then
	cat $name/results/Mutect2/$name.m2.filt.AM.filtered.selected.vcf \
	| java -jar $snpeff_dir/SnpSift.jar filter \
	"( ( FILTER = 'PASS') & (GEN[Tumor].AF >= 0.05) & \
	( ( GEN[Tumor].AD[0] + GEN[Tumor].AD[1]) >= 5 ) & \
	( ( GEN[Normal].AD[0] + GEN[Normal].AD[1]) >= 5 ) & \
	(GEN[Tumor].AD[1] >= 2) & (GEN[Normal].AD[1] <= 1) )" \
	> $name/results/Mutect2/$name.m2.postprocessed.vcf
elif [ $filtering = 'hard' ]; then
	cat $name/results/Mutect2/$name.m2.filt.AM.filtered.selected.vcf \
	| java -jar $snpeff_dir/SnpSift.jar filter \
	"( ( FILTER = 'PASS') & (GEN[Tumor].AF >= 0.1) & \
	( ( GEN[Tumor].AD[0] + GEN[Tumor].AD[1]) >= 10 ) & \
	( ( GEN[Normal].AD[0] + GEN[Normal].AD[1]) >= 10 ) & \
	(GEN[Tumor].AD[1] >= 3) & (GEN[Normal].AD[1] = 0) )" \
	> $name/results/Mutect2/$name.m2.postprocessed.vcf
elif [ $filtering = 'none' ]; then
	cat $name/results/Mutect2/$name.m2.filt.AM.filtered.selected.vcf \
	| java -jar $snpeff_dir/SnpSift.jar filter \
	"( ( FILTER = 'PASS' ) )" \
	> $name/results/Mutect2/$name.m2.postprocessed.vcf
fi

echo '---- Mutect2 Postprocessing II (Filtering out known SNV/Indel using dbSNP or the Sanger Mouse database) ----' | tee -a $name/results/QC/$name.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

bgzip -f $name/results/Mutect2/$name.m2.postprocessed.vcf

tabix -p vcf $name/results/Mutect2/$name.m2.postprocessed.vcf.gz

if [ $filtering = 'all' ]; then
	bcftools isec -C -c none -O z -w 1 \
	-o $name/results/Mutect2/$name.m2.postprocessed.snp_removed.vcf.gz \
	$name/results/Mutect2/$name.m2.postprocessed.vcf.gz \
	$snp_file
elif [ $filtering = 'hard' ]; then
	bcftools isec -C -c none -O z -w 1 \
	-o $name/results/Mutect2/$name.m2.postprocessed.snp_removed.vcf.gz \
	$name/results/Mutect2/$name.m2.postprocessed.vcf.gz \
	$alternate_snp_file
elif [ $filtering = 'none' ]; then
	cp $name/results/Mutect2/$name.m2.postprocessed.vcf.gz $name/results/Mutect2/$name.m2.postprocessed.snp_removed.vcf.gz
fi

bcftools norm -m -any \
$name/results/Mutect2/$name.m2.postprocessed.snp_removed.vcf.gz \
-O z -o $name/results/Mutect2/$name.Mutect2.vcf.gz

gunzip -f $name/results/Mutect2/$name.Mutect2.vcf.gz

echo '---- Mutect2 Postprocessing III (Annotate calls) ----' | tee -a $name/results/QC/$name.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

if [ $species = 'Human' ]; then

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar annotate \
	$dbsnp_file $name/results/Mutect2/$name.Mutect2.vcf \
	> $name/results/Mutect2/$name.Mutect2.ann1.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar \
	annotate $cosmiccoding_file $name/results/Mutect2/$name.Mutect2.ann1.vcf \
	> $name/results/Mutect2/$name.Mutect2.ann2.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar  \
	annotate $cosmicnoncoding_file $name/results/Mutect2/$name.Mutect2.ann2.vcf \
	> $name/results/Mutect2/$name.Mutect2.ann3.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar  \
	annotate $clinvar_file $name/results/Mutect2/$name.Mutect2.ann3.vcf \
	> $name/results/Mutect2/$name.Mutect2.ann4.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar  \
	annotate $gnomadexome_file $name/results/Mutect2/$name.Mutect2.ann4.vcf \
	> $name/results/Mutect2/$name.Mutect2.ann5.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar  \
	annotate $gnomadgenome_file $name/results/Mutect2/$name.Mutect2.ann5.vcf \
	> $name/results/Mutect2/$name.Mutect2.ann6.vcf

	java -Xmx16g -jar "$snpeff_dir"/SnpSift.jar \
	DbNSFP -db $dbnsfp_file -v $name/results/Mutect2/$name.Mutect2.ann6.vcf \
	-f MetaLR_pred,MetaSVM_pred,SIFT_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,PROVEAN_pred \
	> $name/results/Mutect2/"$name".Mutect2.ann7.vcf

	java -Xmx16g -jar $snpeff_dir/snpEff.jar $snpeff_version -canon \
	-csvStats $name/results/Mutect2/$name.Mutect2.annotated.vcf.stats \
	$name/results/Mutect2/"$name".Mutect2.ann7.vcf \
	> $name/results/Mutect2/$name.Mutect2.annotated.vcf

	cat $name/results/Mutect2/$name.Mutect2.annotated.vcf \
	| $snpeff_dir/scripts/vcfEffOnePerLine.pl \
	> $name/results/Mutect2/$name.Mutect2.annotated.one.vcf

	java -jar $snpeff_dir/SnpSift.jar extractFields \
	$name/results/Mutect2/$name.Mutect2.annotated.one.vcf \
	 CHROM POS REF ALT "GEN[Tumor].AF" "GEN[Tumor].AD[0]" "GEN[Tumor].AD[1]" \
	 "GEN[Normal].AD[0]" "GEN[Normal].AD[1]" ANN[*].GENE  ANN[*].EFFECT \
	 ANN[*].IMPACT ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
	 dbNSFP_MetaLR_pred dbNSFP_MetaSVM_pred ID G5 AC AN AF CNT_Coding \
	 CNT_NonCoding CLNDN CLNSIG CLNREVSTAT dbNSFP_SIFT_pred \
	 dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_pred dbNSFP_PROVEAN_pred \
	 > $name/results/Mutect2/"$name".Mutect2.txt

elif [ $species = 'Mouse' ]; then

	java -Xmx16g -jar $snpeff_dir/snpEff.jar $snpeff_version -canon \
	-csvStats $name/results/Mutect2/$name.Mutect2.annotated.vcf.stats \
	$name/results/Mutect2/$name.Mutect2.vcf \
	> $name/results/Mutect2/$name.Mutect2.annotated.vcf

	cat $name/results/Mutect2/$name.Mutect2.annotated.vcf \
	| $snpeff_dir/scripts/vcfEffOnePerLine.pl \
	> $name/results/Mutect2/$name.Mutect2.annotated.one.vcf

	java -jar $snpeff_dir/SnpSift.jar extractFields \
	$name/results/Mutect2/$name.Mutect2.annotated.one.vcf \
	CHROM POS REF ALT "GEN[Tumor].AF" "GEN[Tumor].AD[0]" "GEN[Tumor].AD[1]" \
	"GEN[Normal].AD[0]" "GEN[Normal].AD[1]" ANN[*].GENE  ANN[*].EFFECT \
	ANN[*].IMPACT ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
	> $name/results/Mutect2/$name.Mutect2.txt

fi

sh $repository_dir/SNV_CleanUp.sh $name MS