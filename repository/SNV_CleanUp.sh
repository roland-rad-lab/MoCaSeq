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

	rm $name/results/Mutect2/$name.m2.read-orientation-model.tar.gz
	rm $name/results/Mutect2/$name.m2.postprocessed.snp_removed.vcf.gz
	rm $name/results/Mutect2/$name.m2.filt.selected.vcf.idx
	rm $name/results/Mutect2/$name.m2.filt.selected.vcf
	rm $name/results/Mutect2/$name.m2.filt.vcf.idx
	rm $name/results/Mutect2/$name.m2.filt.vcf
	rm $name/results/Mutect2/$name.m2.postprocessed.vcf.gz
	rm $name/results/Mutect2/$name.m2.postprocessed.vcf.gz.tbi
	rm $name/results/Mutect2/$name.Mutect2.annotated.one.vcf
	rm $name/results/Mutect2/$name.Mutect2.annotated.vcf
	rm $name/results/Mutect2/$name.Mutect2.annotated.vcf.stats
	rm $name/results/Mutect2/$name.Mutect2.annotated.vcf.stats.genes.txt

	if [ $species = 'Human' ]; then
	rm $name/results/Mutect2/$name.Mutect2.ann1.vcf
	rm $name/results/Mutect2/$name.Mutect2.ann2.vcf
	rm $name/results/Mutect2/$name.Mutect2.ann3.vcf
	rm $name/results/Mutect2/$name.Mutect2.ann4.vcf
	rm $name/results/Mutect2/$name.Mutect2.ann5.vcf
	rm $name/results/Mutect2/$name.Mutect2.ann6.vcf
	rm $name/results/Mutect2/$name.Mutect2.ann7.vcf
	fi

elif [ $analysis = 'SS' ]; then

	rm $name/results/Mutect2/$name.$type.m2.postprocessed.snp_removed.vcf.gz
	rm $name/results/Mutect2/$name.$type.Mutect2.annotated.one.vcf
	rm $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf.stats.genes.txt
	rm $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf.stats
	rm $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf
	rm $name/results/Mutect2/$name.$type.m2.postprocessed.vcf.gz
	rm $name/results/Mutect2/$name.$type.m2.postprocessed.vcf.gz.tbi
	rm $name/results/Mutect2/$name.$type.m2.filt.AM.filtered.selected.vcf.idx
	rm $name/results/Mutect2/$name.$type.m2.filt.AM.filtered.selected.vcf
	rm $name/results/Mutect2/$name.$type.m2.filt.AM.vcf.idx
	rm $name/results/Mutect2/$name.$type.m2.filt.AM.vcf
	rm $name/results/Mutect2/$name.$type.m2.filt.AM.vcf.summary
	rm $name/results/Mutect2/$name.$type.m2.filt.AM.filtered.vcf
	rm $name/results/Mutect2/$name.$type.m2.filt.vcf.idx
	rm $name/results/Mutect2/$name.$type.m2.filt.vcf
	rm $name/results/Mutect2/$name.$type.m2.vcf.idx

	if [ $species = 'Human' ]; then
	rm $name/results/Mutect2/$name.$type.Mutect2.ann1.vcf
	rm $name/results/Mutect2/$name.$type.Mutect2.ann2.vcf
	rm $name/results/Mutect2/$name.$type.Mutect2.ann3.vcf
	rm $name/results/Mutect2/$name.$type.Mutect2.ann4.vcf
	rm $name/results/Mutect2/$name.$type.Mutect2.ann5.vcf
	rm $name/results/Mutect2/$name.$type.Mutect2.ann6.vcf
	rm $name/results/Mutect2/$name.$type.Mutect2.ann7.vcf
	fi

elif [ $analysis = 'Manta' ]; then

	rm $name/results/Manta/$name.man.vcf
	rm $name/results/Manta/$name.Manta.annotated.one.vcf
	rm $name/results/Manta/$name.Manta.annotated.vcf
	rm $name/results/Manta/$name.Manta.annotated.vcf.stats
	rm $name/results/Manta/$name.Manta.annotated.vcf.stats.genes.txt
	rm $name/results/Manta/$name.Manta.vcf.gz.tbi

	for type in Tumor Normal;
	do
		rm $name/results/Manta/$name.$type.Manta.annotated.one.vcf
		rm $name/results/Manta/$name.$type.Manta.annotated.vcf
		rm $name/results/Manta/$name.$type.Manta.annotated.vcf.stats
		rm $name/results/Manta/$name.$type.Manta.annotated.vcf.stats.genes.txt
		rm $name/results/Manta/$name.$type.Manta.vcf.gz.tbi
	done

elif [ $analysis = 'Strelka' ]; then
	rm $name/results/Strelka/$name.str.indel.filtered.vcf
	rm $name/results/Strelka/$name.str.indel.postprocessed.vcf.gz
	rm $name/results/Strelka/$name.str.indel.postprocessed.vcf.gz.tbi
	rm $name/results/Strelka/$name.str.indel.postprocessed.vcf.idx
	rm $name/results/Strelka/$name.str.snp.postprocessed.vcf.gz
	rm $name/results/Strelka/$name.str.snp.postprocessed.vcf.gz.tbi
	rm $name/results/Strelka/$name.Strelka.cut
	rm $name/results/Strelka/$name.Strelka.head
	rm $name/results/Strelka/$name.Strelka.tmp
	rm $name/results/Strelka/$name.Strelka.ann1.vcf
	rm $name/results/Strelka/$name.Strelka.ann2.vcf
	rm $name/results/Strelka/$name.Strelka.ann3.vcf
	rm $name/results/Strelka/$name.Strelka.ann4.vcf
	rm $name/results/Strelka/$name.Strelka.ann5.vcf
	rm $name/results/Strelka/$name.Strelka.ann6.vcf
	rm $name/results/Strelka/$name.Strelka.ann7.vcf
	rm $name/results/Strelka/$name.Strelka.annotated.one.vcf
	rm $name/results/Strelka/$name.Strelka.annotated.vcf
	rm $name/results/Strelka/$name.Strelka.annotated.vcf.stats
	rm $name/results/Strelka/$name.Strelka.annotated.vcf.stats.genes.txt
	rm $name/results/Strelka/$name.Strelka.indel.maf.vcf.gz
	rm $name/results/Strelka/$name.Strelka.indel.maf.vcf.gz.tbi
	rm $name/results/Strelka/$name.Strelka.indel.vcf
	rm $name/results/Strelka/$name.Strelka.indel.vep.maf
	rm $name/results/Strelka/$name.Strelka.indel.vep.pairs.tsv
	rm $name/results/Strelka/$name.Strelka.indel.vep.vcf
	rm $name/results/Strelka/$name.Strelka.snp.maf.vcf.gz
	rm $name/results/Strelka/$name.Strelka.snp.maf.vcf.gz.tbi
	rm $name/results/Strelka/$name.Strelka.snp.vcf
	rm $name/results/Strelka/$name.Strelka.snp.vep.maf
	rm $name/results/Strelka/$name.Strelka.snp.vep.pairs.tsv
	rm $name/results/Strelka/$name.Strelka.snp.vep.vcf
	rm $name/results/Strelka/"$name"_vs_NORMAL.vcf
fi

rm -f snpEff_summary.html
rm -f Mutect2FilteringStats.tsv
