

process cnv_kit_matched {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/CNVKit", mode: "copy", pattern: "CNVKit/matched/${meta.sampleName}.matched.CNVKit.*", saveAs: { it.replaceFirst ("^CNVKit/matched/","") }

	input:
		val (genome_build)
		val (reference)
		val (reference_index)
		val (ref_flat)
		tuple path (interval_bed), path (interval_bed_index)
		tuple val (meta), path (bam_normal), path (bai_normal), path (bam_tumor), path (bai_tumor)

	output:
		tuple val (meta), val ("matched"), val ("CNVKit"), env (RESOLUTION), path ("CNVKit/matched/${meta.sampleName}.matched.CNVKit.cnr"), path ("CNVKit/matched/${meta.sampleName}.matched.CNVKit.call.cns"), emit: result
		tuple val (meta), val ("matched"), val ("CNVKit"), env (RESOLUTION), path ("CNVKit/matched/${meta.sampleName}.matched.CNVKit.cns"), emit: cns
		tuple val (meta), val ("matched"), val ("CNVKit"), env (RESOLUTION), path ("CNVKit/matched/${meta.sampleName}.matched.CNVKit.call.cns"), emit: call

	script:
	"""#!/usr/bin/env bash
source ${params.script_base}/file_handling.sh
temp_file_b=\$(moc_mktemp_file .)
trap "rm \${temp_file_b}" EXIT

extract_if_zip ${interval_bed} interval_bed_extracted \${temp_file_b}
mkdir -p CNVKit/matched

touch ${reference_index}
touch ${bai_normal}
touch ${bai_tumor}

cnvkit.py batch \\
	${bam_tumor} \\
	--normal ${bam_normal} \\
	--fasta ${reference} \\
	--output-reference Reference.matched.cnn \\
	--output-dir CNVKit/matched/ \\
	--short-names \\
	--diagram \\
	--scatter \\
	--annotate ${ref_flat} \\
	--access "" \\
	--targets \${interval_bed_extracted} \\
	--drop-low-coverage \\
	-m wgs \\
	-p ${params.cnv_kit.threads}

# The output file path is based on the bam name, lets fix that
# .baseName ~ A nextflow groovy file extension
for file in \$(find CNVKit/matched/ -type f -name "*${bam_normal.baseName}*");
do
		file_new=\$(echo \${file} | sed -e "s/${bam_normal.baseName}/${meta.sampleName}.matched.CNVKit.Normal/")
		echo "Rename '\${file}' to '\${file_new}'"
		mv \${file} \${file_new}
done

for file in \$(find CNVKit/matched/ -type f -name "*${bam_tumor.baseName}*");
do
		file_new=\$(echo \${file} | sed -e "s/${bam_tumor.baseName}/${meta.sampleName}.matched.CNVKit/")
		echo "Rename '\${file}' to '\${file_new}'"
		mv \${file} \${file_new} || true
done

RESOLUTION=\$(cat CNVKit/matched/${meta.sampleName}.matched.CNVKit.targetcoverage.cnn | sed -e '1d;' | awk 'BEGIN{FS=OFS="\\t";total=0;}{total=total+\$3-\$2;}END{print int(total/NR);}')
echo "RESOLUTION: '\${RESOLUTION}'"

	"""

	stub:
	"""#!/usr/bin/env bash
mkdir -p CNVKit/matched

if [[ "${params.stub_json_map?.cnv_kit_matched}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.matched.CNVKit.cnr CNVKit/matched/
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.matched.CNVKit.cns CNVKit/matched/
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.matched.CNVKit.call.cns CNVKit/matched/
fi
touch CNVKit/matched/${meta.sampleName}.matched.CNVKit.cnr
touch CNVKit/matched/${meta.sampleName}.matched.CNVKit.cns
touch CNVKit/matched/${meta.sampleName}.matched.CNVKit.call.cns
RESOLUTION=1000

	"""

}

process cnv_kit_single {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/CNVKit", mode: "copy", pattern: "CNVKit/single/${meta.sampleName}.${type}.CNVKit.*", saveAs: { it.replaceFirst ("^CNVKit/single/","") }

	input:
		val (genome_build)
		val (reference)
		val (reference_index)
		val (ref_flat)
		tuple path (interval_bed), path (interval_bed_index)
		tuple val (meta), val (type), path (bam), path (bai)

	output:
		tuple val (meta), val (type), val ("CNVKit"), env (RESOLUTION), path ("CNVKit/single/${meta.sampleName}.${type}.CNVKit.cnr"), path ("CNVKit/single/${meta.sampleName}.${type}.CNVKit.call.cns"), emit: result
		tuple val (meta), val (type), val ("CNVKit"), env (RESOLUTION), path ("CNVKit/single/${meta.sampleName}.${type}.CNVKit.cns"), emit: cns
		tuple val (meta), val (type), val ("CNVKit"), env (RESOLUTION), path ("CNVKit/single/${meta.sampleName}.${type}.CNVKit.call.cns"), emit: call

	script:
	"""#!/usr/bin/env bash
source ${params.script_base}/file_handling.sh
temp_file_b=\$(moc_mktemp_file .)
trap "rm \${temp_file_b}" EXIT

extract_if_zip ${interval_bed} interval_bed_extracted \${temp_file_b}
mkdir -p CNVKit/single

touch ${reference_index}
touch ${bai}

cnvkit.py batch \\
	${bam} \\
	--normal \\
	--fasta ${reference} \\
	--output-reference Reference.${type}.cnn \\
	--output-dir CNVKit/single/ \\
	--short-names \\
	--diagram \\
	--scatter \\
	--annotate ${ref_flat} \\
	--access "" \\
	--targets \${interval_bed_extracted} \\
	--drop-low-coverage \\
	-m wgs \\
	-p ${params.cnv_kit.threads}

# The output file path is based on the bam name, lets fix that
# .baseName ~ A nextflow groovy file extension
for file in \$(find CNVKit/single/ -type f -name "*${bam.baseName}*");
do
		file_new=\$(echo \${file} | sed -e "s/${bam.baseName}/${meta.sampleName}.${type}.CNVKit/")
		echo "Rename '\${file}' to '\${file_new}'"
		mv \${file} \${file_new} || true
done

RESOLUTION=\$(cat CNVKit/single/${meta.sampleName}.${type}.CNVKit.targetcoverage.cnn | sed -e '1d;' | awk 'BEGIN{FS=OFS="\\t";total=0;}{total=total+\$3-\$2;}END{print int(total/NR);}')
echo "RESOLUTION: '\${RESOLUTION}'"

	"""

	stub:
	"""#!/usr/bin/env bash
mkdir -p CNVKit/single

if [[ "${params.stub_json_map?.cnv_kit_single}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.${type}.CNVKit.cnr CNVKit/single/
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.${type}.CNVKit.cns CNVKit/single/
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.${type}.CNVKit.call.cns CNVKit/single/
fi
touch CNVKit/single/${meta.sampleName}.${type}.CNVKit.cnr
touch CNVKit/single/${meta.sampleName}.${type}.CNVKit.cns
touch CNVKit/single/${meta.sampleName}.${type}.CNVKit.call.cns
RESOLUTION=1000

	"""
}

process cnv_kit_target_bed {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}_PON", mode: "copy"

	input:
		val (genome_build)
		val (reference)
		val (intervals)
		tuple val (meta), val (type), path (bam), path (bai)

	output:
		tuple path ("${genome_build}.target.bed"), path ("${genome_build}.resolution.json"), emit: result	

	script:
	"""#!/usr/bin/env python3.7

from cnvlib import access, autobin, target
from skgenome import tabio

import json

chromosomes = set ("${intervals}".split (","))

access_arr = access.do_access ("${reference}",skip_noncanonical=False).filter (func=lambda x: x["chromosome"] in chromosomes)
tabio.write (access_arr, "access.bed", "bed3")

autobin_args = ['wgs', None, access_arr]
(wgs_depth, target_avg_size), _ = autobin.do_autobin ("${bam}", *autobin_args, bp_per_bin=50000., fasta="${reference}")

print ("wgs_depth: %i\\ntarget_avg_size: %i" % (wgs_depth, target_avg_size))
with open ("${genome_build}.resolution.json", "w") as resolution_file:
	resolution_file.write (json.dumps ({"target_avg_size":target_avg_size,"wgs_depth":wgs_depth}))

annotate = None
short_names = False

target_arr = target.do_target (access_arr, annotate, short_names, True, **({'avg_size': target_avg_size}))
tabio.write(target_arr, "${genome_build}.target.bed", 'bed4')

	"""

	stub:
	"""#!/usr/bin/env bash
if [[ "${params.stub_json_map?.cnv_kit_target_bed}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}_PON/${genome_build}.target.bed .
	cp ${params.stub_dir}/${genome_build}_PON/${genome_build}.resolution.json .
fi

touch ${genome_build}.target.bed
touch ${genome_build}.resolution.json
	"""
}

process cnv_kit_reference {

	publishDir "${params.output_base}/${genome_build}_PON", mode: "copy"

	input:
		val (genome_build)
		val (reference)
		path ("*")
		path (normal_coverage_tsv)

	output:
		path ("${genome_build}.reference.cnn"), emit: result

	script:
	"""#!/usr/bin/env bash

cnvkit.py reference \\
	--fasta ${reference} \\
	--output ${genome_build}.reference.cnn \\
	--no-edge \\
	*.cnn
	"""
}

process cnv_kit_coverage {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/CNVKit", mode: "copy"

	input:
		val (genome_build)
		val (reference)
		tuple path (interval_bed), val (resolution), val (depth)
		tuple val (meta), val (type), path (bam), path (bai)

	output:
		tuple val (meta), val (type), val (resolution), path ("${meta.sampleName}.${type}.coverage.${resolution}.cnn"), emit: result

	script:
	"""#!/usr/bin/env bash
source ${params.script_base}/file_handling.sh
temp_file_b=\$(moc_mktemp_file . bed)
trap "rm \${temp_file_b}" EXIT

extract_if_zip ${interval_bed} interval_bed_extracted \${temp_file_b}

# Stop CNVKit indexing the bam if bam is newer than bai
touch ${bai}

# Giving CNVKit a bed file will give you ratios for those regions
# so one per chromosome if thats what you give
cnvkit.py coverage \\
	--fasta ${reference} \\
	--output ${meta.sampleName}.${type}.coverage.${resolution}.cnn \\
	--processes ${params.cnv_kit.threads} \\
	${bam} \\
	${interval_bed}
	"""

	stub:
	"""#!/usr/bin/env bash
if [[ "${params.stub_json_map?.cnv_kit_coverage}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.${type}.coverage.${resolution}.cnn .
fi
touch ${meta.sampleName}.${type}.coverage.${resolution}.cnn
	"""
}

process cnv_kit_fix {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/CNVKit", mode: "copy"

	input:
		val (genome_build)
		path (reference_cnn)
		tuple val (meta), val (type), val (resolution), path (coverage_cnn) 

	output:
		tuple val (meta), val (type), val("cnv-kit-pon"), val(resolution), path ("${meta.sampleName}.${type}.ratio.${resolution}.cnr"), emit: cnr

	script:
	"""#!/usr/bin/env bash
touch empty.antitargetcoverage.cnn
cnvkit.py fix \\
	--no-edge \\
	--output ${meta.sampleName}.${type}.ratio.${resolution}.cnr \\
	${coverage_cnn} \\
	empty.antitargetcoverage.cnn \\
	${reference_cnn}

	"""
}

process cnv_kit_segment {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/CNVKit", mode: "copy", pattern: "*.cns"

	input:
		val (genome_build)
		each (centre)
		tuple val (meta), val (type), val (coverage_source), val (resolution), path (coverage_cnr)

	output:
		tuple val (meta), val (type), val ("CNVKit.${coverage_source}.${centre}"), val (resolution), path ("${meta.sampleName}.${type}.CNVKit.${coverage_source}.${centre}.cnr"), path ("${meta.sampleName}.${type}.CNVKit.${coverage_source}.${centre}.call.cns"), emit: result
		tuple val (meta), val (type), val ("CNVKit.${coverage_source}.${centre}"), val (resolution), path ("${meta.sampleName}.${type}.CNVKit.${coverage_source}.${centre}.cnr"), emit: cnr
		tuple val (meta), val (type), val ("CNVKit.${coverage_source}.${centre}"), val (resolution), path ("${meta.sampleName}.${type}.CNVKit.${coverage_source}.${centre}.call.cns"), emit: call
		path ("${meta.sampleName}.${type}.CNVKit.${coverage_source}.${centre}.cns")

	script:
	"""#!/usr/bin/env bash

cnvkit.py segment \\
	-p ${params.cnv_kit.threads} \\
	--drop-low-coverage \\
	-o ${meta.sampleName}.${type}.CNVKit.${coverage_source}.${centre}.cns \\
	${coverage_cnr}

cnvkit.py call \\
	-m none \\
	--center ${centre} \\
	-o ${meta.sampleName}.${type}.CNVKit.${coverage_source}.${centre}.call.cns \\
	${meta.sampleName}.${type}.CNVKit.${coverage_source}.${centre}.cns

cnvkit.py call \\
	-m none \\
	--center ${centre} \\
	-o ${meta.sampleName}.${type}.CNVKit.${coverage_source}.${centre}.cnr
	${coverage_cnr}

	"""

	stub:
	"""#!/usr/bin/env bash
if [[ "${params.stub_json_map?.cnv_kit_segment}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.${type}.CNVKit.${coverage_source}.${centre}.cnr .
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.${type}.CNVKit.${coverage_source}.${centre}.cns .
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.${type}.CNVKit.${coverage_source}.${centre}.call.cns .
fi

touch ${meta.sampleName}.${type}.CNVKit.${coverage_source}.${centre}.cnr
touch ${meta.sampleName}.${type}.CNVKit.${coverage_source}.${centre}.cns
touch ${meta.sampleName}.${type}.CNVKit.${coverage_source}.${centre}.call.cns
	"""
}

process cnv_kit_plot {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/CNVKit", mode: "copy"

	input:
		val (genome_build)
		tuple path (interval_bed), path (interval_bed_index)
		tuple val (meta), val (type), val (coverage_source), val (resolution), path (cnr_file), path (call_file)

	output:
		path ("${meta.sampleName}.*.pdf")

	script:
	"""#!/usr/bin/env Rscript

library (dplyr)
library (ggplot2)
library (gridExtra)

interval_file <- gzfile ("${interval_bed}", 'rt')
data_interval <- read.table (file=interval_file,sep="\\t",header=F,stringsAsFactors=F)
names (data_interval) <- c("Chrom", "Start", "End")
head (data_interval)

data_ratio <- read.table (file="${cnr_file}",sep="\\t",header=T,stringsAsFactors=F)
head (data_ratio)

data_call <- read.table (file="${call_file}",sep="\\t",header=T,stringsAsFactors=F)
head (data_call)

data_interval_plot <- data_interval %>%
	dplyr::mutate (Chrom=as.character (Chrom),End=as.numeric (End)) %>%
	dplyr::mutate (CumulativeStart=cumsum (End)-End) %>%
	dplyr::mutate (CumulativeEnd=cumsum (End)) %>%
	dplyr::mutate (CumulativeMidpoint=(CumulativeStart+CumulativeEnd)/2) %>%
	data.frame

data_ratio_plot <- data_ratio %>%
	dplyr::filter (!is.na (log2)) %>%
	dplyr::mutate (midpoint=(start+end)/2) %>%
	dplyr::inner_join (data_interval_plot,by=c("chromosome"="Chrom"),suffix=c("",".Chrom")) %>%
	dplyr::mutate (Start.Genome=start+CumulativeStart,End.Genome=end+CumulativeStart,Midpoint.Genome=midpoint+CumulativeStart) %>%
	data.frame

data_call_plot <- data_call %>%
	dplyr::inner_join (data_interval_plot,by=c("chromosome"="Chrom"),suffix=c("",".Chrom")) %>%
	dplyr::mutate (Start.Genome=start+CumulativeStart,End.Genome=end+CumulativeStart) %>%
	data.frame

head (data_interval_plot)
head (data_ratio_plot)
head (data_call_plot)

pdf (file="${meta.sampleName}.${type}.${coverage_source}.${resolution}.genome.pdf",width=16,height=4)

ggplot (data_ratio_plot) +
	#geom_segment (aes(x=Start.Genome,y=log2,xend=End.Genome,yend=log2),colour="#636363") +
	geom_point (aes(x=Midpoint.Genome,y=log2),shape=".",colour="#c7c7c7") +
	geom_segment (data=data_call_plot,aes(x=Start.Genome,y=log2,xend=End.Genome,yend=log2),colour="red") +
	geom_vline (data=data_interval_plot,aes(xintercept=CumulativeStart),colour="#D3D3D3") +
	geom_text (data=data_interval_plot,aes(x=CumulativeMidpoint,y=2.1,label=Chrom),size=2) +
	coord_cartesian (xlim=c(0,data_interval_plot %>% pull (CumulativeEnd) %>% max ()), ylim=c(-2,2),expand=F,clip="off") +
	xlab ("Genome") +
	theme_bw () +
	theme (
		panel.grid.major.x=element_blank (),
		panel.grid.minor.x=element_blank (),
		axis.ticks.x=element_blank (),
		axis.text.x=element_blank (),
		plot.margin = unit(c(1,0.5,0.5,0.5), "cm")
	)

dev.off ()

chromosomes <- data_interval %>% pull (Chrom)
plot_list <- vector ("list",length (chromosomes))


for ( i in seq_along (chromosomes) )
{
	plot_list[[i]] <- ggplot (data_ratio_plot %>% filter (chromosome==!!chromosomes[[i]]) %>% data.frame) +
		#geom_segment (aes(x=start,y=log2,xend=end,yend=log2),colour="#636363") +
		geom_point (aes(x=midpoint,y=log2),shape=".",colour="#c7c7c7") +
		geom_segment (data=data_call_plot %>% filter (chromosome==!!chromosomes[[i]]) %>% data.frame,aes(x=start,y=log2,xend=end,yend=log2),colour="red") +
		scale_x_continuous (labels=scales::number_format (big.mark=",",scale=1e-06,suffix=" Mb",accuracy=0.1)) +
		coord_cartesian (xlim=c(0,data_interval_plot %>% filter (Chrom==!!chromosomes[[i]]) %>% pull (End)),ylim=c(-2,2)) +
		labs (title=chromosomes[[i]]) +
		theme_bw () +
		theme (
			panel.grid.major.x=element_blank (),
			panel.grid.minor.x=element_blank (),
			plot.margin=unit (c(5.5,25.5,5.5,5.5),"pt")
		)
}

# Need to do this outside of pdf call to prevent blank first page
p <- marrangeGrob (plot_list,nrow=1,ncol=1)

pdf (file="${meta.sampleName}.${type}.${coverage_source}.${resolution}.chromosomes.pdf",width=9)
p
dev.off ()


	"""

	stub:
	"""#!/usr/bin/env bash

touch ${meta.sampleName}.${type}.${coverage_source}.${resolution}.chromosomes.pdf
touch ${meta.sampleName}.${type}.${coverage_source}.${resolution}.genome.pdf

	"""
}

