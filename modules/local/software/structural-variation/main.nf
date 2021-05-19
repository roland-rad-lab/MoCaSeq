



process structural_variation_matched {
	tag "${meta.sampleName}"

	input:
		tuple path (interval_bed), path(interval_bed_index)
		tuple val (meta), path (manta_vcf), path (manta_vcf_index), path (delly_bcf), path (normal_cns), path (tumor_cns)

	output:
		tuple val (meta), path ("${meta.sampleName}.manta.tsv"), path ("${meta.sampleName}.delly.tsv"), emit: result

	script:
	"""#!/usr/bin/env bash

echo -e "SAMPLE\\tCHROM\\tPOS\\tID\\tREF\\tALT\\tINFO.IMPRECISE\\tINFO.SVTYPE\\tINFO.SVLEN\\tINFO.END\\tINFO.CIPOS\\tINFO.CIEND\\tINFO.CIGAR\\tINFO.MATEID\\tINFO.EVENT\\tINFO.HOMLEN\\tINFO.HOMSEQ\\tINFO.SVINSLEN\\tINFO.SVINSSEQ\\tINFO.LEFT_SVINSSEQ\\tINFO.RIGHT_SVINSSEQ\\tINFO.BND_DEPTH\\tINFO.MATE_BND_DEPTH\\tINFO.SOMATIC\\tINFO.SOMATICSCORE\\tINFO.JUNCTION_SOMATICSCORE\\tPR\\tSR" | tr "[A-Z]" "[a-z]" > ${meta.sampleName}.manta.tsv
bcftools view -f PASS ${manta_vcf} | bcftools query -f '[%SAMPLE\\t%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/IMPRECISE\\t%INFO/SVTYPE\\t%INFO/SVLEN\\t%INFO/END\\t%INFO/CIPOS\\t%INFO/CIEND\\t%INFO/CIGAR\\t%INFO/MATEID\\t%INFO/EVENT\\t%INFO/HOMLEN\\t%INFO/HOMSEQ\\t%INFO/SVINSLEN\\t%INFO/SVINSSEQ\\t%INFO/LEFT_SVINSSEQ\\t%INFO/RIGHT_SVINSSEQ\\t%INFO/BND_DEPTH\\t%INFO/MATE_BND_DEPTH\\t%INFO/SOMATIC\\t%INFO/SOMATICSCORE\\t%INFO/JUNCTION_SOMATICSCORE\\t%PR\\t%SR\n]' >> ${meta.sampleName}.manta.tsv

echo -e "SAMPLE\\tCHROM\\tPOS\\tID\\tREF\\tALT\\tINFO.IMPRECISE\\tINFO.SVTYPE\\tINFO.SVLEN\\tINFO.END\\tINFO.CIPOS\\tINFO.CIEND\\tINFO.CHR2\\tINFO.POS2\\tRR\\tRV" | tr "[A-Z]" "[a-z]" > ${meta.sampleName}.delly.tsv
bcftools view -f PASS ${delly_bcf} | bcftools query -f '[%SAMPLE\\t%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/IMPRECISE\\t%INFO/SVTYPE\\t%INFO/SVLEN\\t%INFO/END\\t%INFO/CIPOS\\t%INFO/CIEND\\t%INFO/CHR2\\t%INFO/POS2\\t%RR\\t%RV\n]' >> ${meta.sampleName}.delly.tsv

echo -e "sample\tchrom\tstart\tend\tlog2" > ${meta.sampleName}.cnvkit.tsv
cat ${normal_cns} | sed -e '1d;' | awk BEGIN'{FS=OFS="\\t";}{print "Normal",\$1,\$2,\$3,\$5;}' >> ${meta.sampleName}.cnvkit.tsv
cat ${tumor_cns} | sed -e '1d;' | awk BEGIN'{FS=OFS="\\t";}{print "Tumor",\$1,\$2,\$3,\$5;}' >> ${meta.sampleName}.cnvkit.tsv

	"""
}

process structural_variation_matched_merge {
	tag "${meta.sampleName}"

	input:
		tuple val (meta), path (manta_tsv), path (delly_tsv)


	script:
	"""#!/usr/bin/env Rscript
library (dplyr)
library (stringr)

source ("${params.script_base}/SV_library.R")

manta_tsv_file_path <- "${manta_tsv}"
delly_tsv_file_path <- "${delly_tsv}"

data_manta <- read.table (file=manta_tsv_file_path,header=T,sep="\t",stringsAsFactors=F)
head (data_manta)

data_delly <- read.table (file=delly_tsv_file_path,header=T,sep="\\t",stringsAsFactors=F)
head (data_delly)

data_manta_sv <- sv_from_manta (data_manta)
head (data_manta_sv)

data_delly_sv <- sv_from_delly (data_delly)
head (data_delly_sv)

data_sv <- rbind (data_manta_sv %>% mutate(caller="manta") %>% data.frame,data_delly_sv %>% mutate(caller="delly") %>% data.frame)

write.table (data_sv,file="sv.combined.tsv",sep="\\t",row.names=F,quote=F)

	"""
}


