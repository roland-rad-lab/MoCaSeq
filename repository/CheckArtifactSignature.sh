# NOT WORKING AND UNFINISHED SCRIPT


sudo docker run \
-it --entrypoint=/bin/bash \
-v /home/rad/Downloads/maf_filtering:/var/pipeline/ \
-v /mnt/3TBVol/MoCaSeq_WD/MoCaSeq_ref:/var/pipeline/ref \
-v /home/rad/packages/MoCaSeq:/opt/MoCaSeq \
mocaseq-human

sudo docker run \
-it --entrypoint=/bin/bash \
-v /media/rad/HDD1/niklas_tmp:/var/pipeline/ \
-v /media/rad/SSD1/MoCaSeq_ref:/var/pipeline/ref \
rolandradlab/mocaseq

GATK=4.1.7.0
GATK_dir=$(echo /opt/gatk-${GATK})
genomes_dir=/var/pipeline/ref/
genome_dir=$genomes_dir/GRCh38.p12
genome_file=$genome_dir/GRCh38.p12.fna


# for single sample
name=1N7T85
type=Tumor



java -jar $GATK_dir/gatk.jar FilterMutectCalls \
--variant $name/results/Mutect2/$name.$type.Mutect2.vcf \
--output $name/results/Mutect2/$name.$type.m2.filt.vcf \
--reference $genome_file

if [ $artefact_type = 'none' ]; then
	cp $name/results/Mutect2/$name.$type.m2.filt.vcf \
	$name/results/Mutect2/$name.$type.m2.filt.AM.vcf
elif [ $artefact_type = 'CT' ]; then
	java -jar $GATK_dir/gatk.jar FilterByOrientationBias \
	-V $name/results/Mutect2/$name.$type.m2.filt.vcf -P \
	$name/results/QC/$name.$type.bam.artifacts.pre_adapter_detail_metrics \
	--artifact-modes C/T --output $name/results/Mutect2/$name.$type.m2.filt.AM.vcf
elif [ $artefact_type = 'GT' ]; then
	java -jar $GATK_dir/gatk.jar FilterByOrientationBias \
	-V $name/results/Mutect2/$name.$type.m2.filt.vcf -P \
	$name/results/QC/$name.$type.bam.artifacts.pre_adapter_detail_metrics \
	--artifact-modes G/T --output $name/results/Mutect2/$name.$type.m2.filt.AM.vcf
fi

if [ $artefact_type = 'none' ]; then
	cp $name/results/Mutect2/$name.$type.m2.filt.AM.vcf \
	$name/results/Mutect2/$name.$type.m2.filt.AM.filtered.vcf
elif [ $artefact_type = 'CT' ] || [ $artefact_type = 'GT' ]; then
	cat $name/results/Mutect2/$name.$type.m2.filt.AM.vcf \
	| java -jar $snpeff_dir/SnpSift.jar filter \
	"( ( ( FILTER = 'PASS'  ) & (exists GEN[$type].OBP) & \
	(GEN[$type].OBP <= 0.05) ) | ( ( FILTER = 'PASS' ) ) )" \
	> $name/results/Mutect2/$name.$type.m2.filt.AM.filtered.vcf
fi
