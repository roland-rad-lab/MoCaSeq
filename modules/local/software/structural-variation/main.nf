



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

echo -e "SAMPLE\\tCHROM\\tPOS\\tID\\tREF\\tALT\\tINFO.IMPRECISE\\tINFO.SVTYPE\\tINFO.SVLEN\\tINFO.END\\tINFO.CIPOS\\tINFO.CIEND\\tINFO.CHR2\\tINFO.POS2\\tRR\\tRV\\tGT\\tGQ\\tFT" | tr "[A-Z]" "[a-z]" > ${meta.sampleName}.delly.tsv
bcftools view -f PASS ${delly_bcf} | bcftools query -f '[%SAMPLE\\t%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/IMPRECISE\\t%INFO/SVTYPE\\t%INFO/SVLEN\\t%INFO/END\\t%INFO/CIPOS\\t%INFO/CIEND\\t%INFO/CHR2\\t%INFO/POS2\\t%RR\\t%RV\\t%GT\\t%GQ\\t%FT\n]' >> ${meta.sampleName}.delly.tsv

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
library (IRanges)

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

print ("data_sv")
print (head(data_sv))
print (class(data_sv[,"pos1"]))

# https://www.bioconductor.org/packages/release/bioc/html/StructuralVariantAnnotation.html
# ^ Can read from the VCFs and do the merging, plus output to bedpe and circos plots
# at some point borrow more of their merge ideas (i.e. CI)
# https://github.com/PapenfussLab/StructuralVariantAnnotation

data_sv_merged <- data_sv %>%
	#filter (sample=="Tumor",SVtype=="t2tINV",chrom1=="3",chrom2=="1") %>%
	#filter (sample=="Normal",SVtype=="DUP",chrom1=="14",chrom2=="14") %>%
	#bind_cols (event_id_sample_sv_chrom=group_indices (.,sample,SVtype,chrom1,chrom2)) %>%
	#mutate (event_id_sample_sv_chrom=group_by (.,sample,SVtype,chrom1,chrom2) %>% group_indices ()) %>%
	group_by (sample,SVtype,chrom1,chrom2) %>%
	mutate(event_id_sample_sv_chrom=group()) %>%
	group_modify (sv_find_event) %>%
	ungroup () %>%
	dplyr::select (-c(region1,region2)) %>%
	dplyr::rename (event_id_sample_sv_chrom_event=event_id) %>%
	bind_cols (event_id_sample=group_indices (.,event_id_sample_sv_chrom,event_id_sample_sv_chrom_event)) %>%
	bind_cols (event_id_sv_chrom=group_indices (.,SVtype,chrom1,chrom2)) %>%
	group_by (SVtype,chrom1,chrom2) %>%
	group_modify (sv_find_event) %>%
	ungroup () %>%
	dplyr::select (-c(region1,region2)) %>%
	dplyr::rename (event_id_sv_chrom_event=event_id) %>%
	bind_cols (event_id=group_indices (.,event_id_sv_chrom,event_id_sv_chrom_event)) %>%
	group_by (event_id) %>%
	#group_modify(function (data,group) { print (data[,"sample"]);print ("Normal" %in% data[,"sample"]); cbind (data,somatic="Normal" %in% data[,"sample"]) }) %>%
	mutate (somatic=match("Normal",sample,nomatch=0)) %>%
	ungroup () %>%
	dplyr::select (-c (event_id_sample_sv_chrom,event_id_sample_sv_chrom_event,event_id_sv_chrom,event_id_sv_chrom_event)) %>%
	data.frame

write.table (data_sv_merged,file="${meta.sampleName}.sv.events.tsv",sep="\\t",row.names=F,quote=F)

	"""
}


