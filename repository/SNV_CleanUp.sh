#!/bin/bash

##########################################################################################
##
## SNV_CleanUp.sh
##
## Removes all intermediary files from SNV/Indel and LOH detection.
##
##########################################################################################

name=$1
analysis=$2
type=$3
species=$4

if [ $analysis = 'MS' ]; then

	rm -f $name/results/Mutect2/$name.m2.read-orientation-model.tar.gz
	rm -f $name/results/Mutect2/$name.m2.postprocessed.snp_removed.vcf.gz
	rm -f $name/results/Mutect2/$name.m2.filt.selected.vcf.idx
	rm -f $name/results/Mutect2/$name.m2.filt.selected.vcf
	rm -f $name/results/Mutect2/$name.m2.filt.vcf.idx
	rm -f $name/results/Mutect2/$name.m2.filt.vcf
	rm -f $name/results/Mutect2/$name.m2.postprocessed.vcf.gz
	rm -f $name/results/Mutect2/$name.m2.postprocessed.vcf.gz.tbi
	#rm -f $name/results/Mutect2/$name.Mutect2.annotated.one.vcf
	rm -f $name/results/Mutect2/$name.Mutect2.annotated.vcf
	rm -f $name/results/Mutect2/$name.Mutect2.annotated.vcf.stats
	rm -f $name/results/Mutect2/$name.Mutect2.annotated.vcf.stats.genes.txt

	if [ $species = 'Human' ]; then
	rm -f $name/results/Mutect2/$name.Mutect2.ann1.vcf
	rm -f $name/results/Mutect2/$name.Mutect2.ann2.vcf
	rm -f $name/results/Mutect2/$name.Mutect2.ann3.vcf
	rm -f $name/results/Mutect2/$name.Mutect2.ann4.vcf
	rm -f $name/results/Mutect2/$name.Mutect2.ann5.vcf
	rm -f $name/results/Mutect2/$name.Mutect2.ann6.vcf
	rm -f $name/results/Mutect2/$name.Mutect2.ann7.vcf
	fi

elif [ $analysis = 'SS' ]; then

	rm -f $name/results/Mutect2/$name.$type.m2.postprocessed.snp_removed.vcf.gz
	#rm -f $name/results/Mutect2/$name.$type.Mutect2.annotated.one.vcf
	rm -f $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf.stats.genes.txt
	rm -f $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf.stats
	rm -f $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf
	rm -f $name/results/Mutect2/$name.$type.m2.postprocessed.vcf.gz
	rm -f $name/results/Mutect2/$name.$type.m2.postprocessed.vcf.gz.tbi
	rm -f $name/results/Mutect2/$name.$type.m2.filt.AM.filtered.selected.vcf.idx
	rm -f $name/results/Mutect2/$name.$type.m2.filt.AM.filtered.selected.vcf
	rm -f $name/results/Mutect2/$name.$type.m2.filt.AM.vcf.idx
	rm -f $name/results/Mutect2/$name.$type.m2.filt.AM.vcf
	rm -f $name/results/Mutect2/$name.$type.m2.filt.AM.vcf.summary
	rm -f $name/results/Mutect2/$name.$type.m2.filt.AM.filtered.vcf
	rm -f $name/results/Mutect2/$name.$type.m2.filt.vcf.idx
	rm -f $name/results/Mutect2/$name.$type.m2.filt.vcf
	rm -f $name/results/Mutect2/$name.$type.m2.vcf.idx

	if [ $species = 'Human' ]; then
	rm -f $name/results/Mutect2/$name.$type.Mutect2.ann1.vcf
	rm -f $name/results/Mutect2/$name.$type.Mutect2.ann2.vcf
	rm -f $name/results/Mutect2/$name.$type.Mutect2.ann3.vcf
	rm -f $name/results/Mutect2/$name.$type.Mutect2.ann4.vcf
	rm -f $name/results/Mutect2/$name.$type.Mutect2.ann5.vcf
	rm -f $name/results/Mutect2/$name.$type.Mutect2.ann6.vcf
	rm -f $name/results/Mutect2/$name.$type.Mutect2.ann7.vcf
	fi

elif [ $analysis = 'Manta' ]; then

	rm -f $name/results/Manta/$name.man.vcf
	rm -f $name/results/Manta/$name.Manta.annotated.one.vcf
	rm -f $name/results/Manta/$name.Manta.annotated.vcf
	rm -f $name/results/Manta/$name.Manta.annotated.vcf.stats
	rm -f $name/results/Manta/$name.Manta.annotated.vcf.stats.genes.txt
	rm -f $name/results/Manta/$name.Manta.vcf.gz.tbi

	for type in Tumor Normal;
	do
		rm -f $name/results/Manta/$name.$type.Manta.annotated.one.vcf
		rm -f $name/results/Manta/$name.$type.Manta.annotated.vcf
		rm -f $name/results/Manta/$name.$type.Manta.annotated.vcf.stats
		rm -f $name/results/Manta/$name.$type.Manta.annotated.vcf.stats.genes.txt
		rm -f $name/results/Manta/$name.$type.Manta.vcf.gz.tbi
	done

elif [ $analysis = 'Strelka' ]; then
	rm -f $name/results/Strelka/$name.str.indel.filtered.vcf
	rm -f $name/results/Strelka/$name.str.indel.postprocessed.vcf.gz
	rm -f $name/results/Strelka/$name.str.indel.postprocessed.vcf.gz.tbi
	rm -f $name/results/Strelka/$name.str.indel.postprocessed.vcf.idx
	rm -f $name/results/Strelka/$name.str.snp.postprocessed.vcf.gz
	rm -f $name/results/Strelka/$name.str.snp.postprocessed.vcf.gz.tbi
	rm -f $name/results/Strelka/$name.Strelka.cut
	rm -f $name/results/Strelka/$name.Strelka.head
	rm -f $name/results/Strelka/$name.Strelka.tmp
	rm -f $name/results/Strelka/$name.Strelka.ann1.vcf
	rm -f $name/results/Strelka/$name.Strelka.ann2.vcf
	rm -f $name/results/Strelka/$name.Strelka.ann3.vcf
	rm -f $name/results/Strelka/$name.Strelka.ann4.vcf
	rm -f $name/results/Strelka/$name.Strelka.ann5.vcf
	rm -f $name/results/Strelka/$name.Strelka.ann6.vcf
	rm -f $name/results/Strelka/$name.Strelka.ann7.vcf
	#rm -f $name/results/Strelka/$name.Strelka.annotated.one.vcf
	rm -f $name/results/Strelka/$name.Strelka.annotated.vcf
	rm -f $name/results/Strelka/$name.Strelka.annotated.vcf.stats
	rm -f $name/results/Strelka/$name.Strelka.annotated.vcf.stats.genes.txt
	rm -f $name/results/Strelka/$name.Strelka.indel.maf.vcf.gz
	rm -f $name/results/Strelka/$name.Strelka.indel.maf.vcf.gz.tbi
	rm -f $name/results/Strelka/$name.Strelka.indel.vcf
	rm -f $name/results/Strelka/$name.Strelka.indel.vep.maf
	rm -f $name/results/Strelka/$name.Strelka.indel.vep.pairs.tsv
	rm -f $name/results/Strelka/$name.Strelka.indel.vep.vcf
	rm -f $name/results/Strelka/$name.Strelka.snp.maf.vcf.gz
	rm -f $name/results/Strelka/$name.Strelka.snp.maf.vcf.gz.tbi
	rm -f $name/results/Strelka/$name.Strelka.snp.vcf
	rm -f $name/results/Strelka/$name.Strelka.snp.vep.maf
	rm -f $name/results/Strelka/$name.Strelka.snp.vep.pairs.tsv
	rm -f $name/results/Strelka/$name.Strelka.snp.vep.vcf
	rm -f $name/results/Strelka/"$name"_vs_NORMAL.vcf
fi

rm -f -f snpEff_summary.html
rm -f -f Mutect2FilteringStats.tsv
