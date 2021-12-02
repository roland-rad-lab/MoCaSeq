

process cnv_kit_matched {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/CNVKit", mode: "copy", pattern: "CNVKit/matched/${meta.sampleName}.matched.cns", saveAs: { it.replaceFirst ("^CNVKit/matched/","") }

	input:
		val (genome_build)
		val (reference)
		val (reference_index)
		val (ref_flat)
		tuple path (interval_bed), path (interval_bed_index)
		tuple val (meta), path (bam_normal), path (bai_normal), path (bam_tumor), path (bai_tumor)

	output:
		tuple val (meta), val ("matched"), path ("CNVKit/matched/${meta.sampleName}.matched.cns"), emit: cns

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
		file_new=\$(echo \${file} | sed -e "s/${bam_normal.baseName}/${meta.sampleName}.matched.Normal/")
		echo "Rename '\${file}' to '\${file_new}'"
		mv \${file} \${file_new}
done

for file in \$(find CNVKit/matched/ -type f -name "*${bam_tumor.baseName}*");
do
		file_new=\$(echo \${file} | sed -e "s/${bam_tumor.baseName}/${meta.sampleName}.matched/")
		echo "Rename '\${file}' to '\${file_new}'"
		mv \${file} \${file_new} || true
done

	"""

	stub:
	"""#!/usr/bin/env bash
mkdir -p CNVKit/matched

if [[ "${params.stub_json_map?.cnv_kit_matched}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.matched.cns CNVKit/matched/
fi
touch CNVKit/matched/${meta.sampleName}.matched.cns
	"""

}

process cnv_kit_single {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/CNVKit", mode: "copy", pattern: "CNVKit/single/${meta.sampleName}.${type}.cns", saveAs: { it.replaceFirst ("^CNVKit/single/","") }

	input:
		val (genome_build)
		val (reference)
		val (reference_index)
		val (ref_flat)
		tuple path (interval_bed), path (interval_bed_index)
		tuple val (meta), val (type), path (bam), path (bai)

	output:
		tuple val (meta), val (type), path ("CNVKit/single/${meta.sampleName}.${type}.cns"), emit: cns

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
		file_new=\$(echo \${file} | sed -e "s/${bam.baseName}/${meta.sampleName}.${type}/")
		echo "Rename '\${file}' to '\${file_new}'"
		mv \${file} \${file_new} || true
done

	"""

	stub:
	"""#!/usr/bin/env bash
mkdir -p CNVKit/single

if [[ "${params.stub_json_map?.cnv_kit_single}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.${type}.cns CNVKit/single/
fi
touch CNVKit/single/${meta.sampleName}.${type}.cns
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
		tuple val (meta), val (type), val (coverage_source), val (resolution), path (coverage_cnr)

	output:
		tuple val (meta), val (type), val ("${coverage_source}-cnv-kit"), val (resolution), path ("${meta.sampleName}.${type}.${coverage_source}.cns"), emit: cns
		tuple val (meta), val (type), val ("${coverage_source}-cnv-kit"), val (resolution), path ("${meta.sampleName}.${type}.${coverage_source}.mode.call.cns"), emit: call

	script:
	"""#!/usr/bin/env bash

cnvkit.py segment -o ${meta.sampleName}.${type}.${coverage_source}.cns ${coverage_cnr}
cnvkit.py call -o ${meta.sampleName}.${type}.${coverage_source}.mode.call.cns --center mode ${meta.sampleName}.${type}.${coverage_source}.cns
	"""

	stub:
	"""#!/usr/bin/env bash
if [[ "${params.stub_json_map?.cnv_kit_segment}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.${type}.${coverage_source}.cns .
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.${type}.${coverage_source}.mode.call.cns .
fi

touch ${meta.sampleName}.${type}.${coverage_source}.cns
touch ${meta.sampleName}.${type}.${coverage_source}.mode.call.cns
	"""
}

