#!/bin/bash

##########################################################################################
##
## SNV_MantaPostprocessing.sh
##
## Postprocessing for Manta.
##
##########################################################################################

name=$1
species=$2
config_file=$3
runmode=$4
type=$5

. $config_file

if [ $runmode = 'MS' ]; then

	cp $name/results/Manta/results/variants/somaticSV.vcf.gz $name/results/Manta/$name.man.vcf.gz

	gunzip -f $name/results/Manta/$name.man.vcf.gz

	cat $name/results/Manta/$name.man.vcf | java -jar "$snpeff_dir"/SnpSift.jar filter "( ( ( FILTER = 'PASS' ) & ( ( GEN[Tumor].SR[1] + GEN[Tumor].SR[0] ) * 0.05 <= GEN[Tumor].SR[1] ) ) | ( ( FILTER = 'PASS' ) & ( ( GEN[Tumor].PR[1] + GEN[Tumor].PR[0] ) * 0.05 <= GEN[Tumor].PR[1] ) ) )" > $name/results/Manta/$name.Manta.vcf

	bgzip $name/results/Manta/$name.Manta.vcf

	tabix -p vcf $name/results/Manta/$name.Manta.vcf.gz

	bcftools stats $name/results/Manta/$name.Manta.vcf.gz > $name/results/Manta/$name.Manta.vcf.gz.stats

	gunzip -f $name/results/Manta/$name.Manta.vcf.gz

	java -Xmx16g -jar $snpeff_dir/snpEff.jar $snpeff_version -canon \
	-csvStats $name/results/Manta/$name.Manta.annotated.vcf.stats \
	$name/results/Manta/$name.Manta.vcf \
	> $name/results/Manta/$name.Manta.annotated.vcf

	cat $name/results/Manta/$name.Manta.annotated.vcf | "$snpeff_dir"/scripts/vcfEffOnePerLine.pl > $name/results/Manta/$name.Manta.annotated.one.vcf

	java -jar "$snpeff_dir"/SnpSift.jar extractFields $name/results/Manta/$name.Manta.annotated.one.vcf CHROM POS REF ALT "GEN[Tumor].SR[0]" "GEN[Tumor].SR[1]" "GEN[Tumor].PR[0]" "GEN[Tumor].PR[1]" "GEN[Normal].SR[0]" "GEN[Normal].SR[1]" "GEN[Normal].PR[0]" "GEN[Normal].PR[1]" ANN[*].GENE  ANN[*].EFFECT ANN[*].IMPACT ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P > $name/results/Manta/$name.Manta.txt

elif [ $runmode = 'SS' ]; then

	cp $name/results/Manta/results/variants/diploidSV.vcf.gz $name/results/Manta/$name.man.vcf.gz

	gunzip -f $name/results/Manta/$name.man.vcf.gz

	cat $name/results/Manta/$name.man.vcf | java -jar "$snpeff_dir"/SnpSift.jar filter "( ( FILTER = 'PASS' ) )" > $name/results/Manta/$name.$type.Manta.vcf

	bgzip $name/results/Manta/$name.$type.Manta.vcf

	tabix -p vcf $name/results/Manta/$name.$type.Manta.vcf.gz

	bcftools stats $name/results/Manta/$name.$type.Manta.vcf.gz > $name/results/Manta/$name.$type.Manta.vcf.gz.stats

	gunzip -f $name/results/Manta/$name.$type.Manta.vcf.gz

	java -Xmx16g -jar "$snpeff_dir"/snpEff.jar $snpeff_version -canon \
	-csvStats $name/results/Manta/$name.Manta.annotated.vcf.stats \
	$name/results/Manta/$name.$type.Manta.vcf \
	> $name/results/Manta/$name.$type.Manta.annotated.vcf

	cat $name/results/Manta/$name.$type.Manta.annotated.vcf | "$snpeff_dir"/scripts/vcfEffOnePerLine.pl > $name/results/Manta/$name.$type.Manta.annotated.one.vcf

	java -jar "$snpeff_dir"/SnpSift.jar extractFields $name/results/Manta/$name.$type.Manta.annotated.one.vcf CHROM POS REF ALT "GEN[$type].SR[0]" "GEN[$type].SR[1]" "GEN[$type].PR[0]" "GEN[$type].PR[1]" ANN[*].GENE  ANN[*].EFFECT ANN[*].IMPACT ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P > $name/results/Manta/$name.$type.Manta.txt
fi

# python2 $svtools_dir/vcfToBedpe \
# -i $name/results/Manta/$name.$type.Manta.vcf \
# -o $name/results/Manta/$name.$type.Manta.bedpe

sh $repository_dir/SNV_CleanUp.sh $name Manta
