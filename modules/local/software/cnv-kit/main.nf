

process cnv_kit_matched {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/CNVKit", mode: "copy", pattern: "CNVKit/matched/${meta.sampleName}.matched.cns", saveAs: { it.replaceFirst ("^CNVKit/matched/","") }

	input:
		val (genome_build)
		val (reference)
		val (reference_flat)
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
cnvkit.py batch \\
	${bam_tumor} \\
	--normal ${bam_normal} \\
	--fasta ${reference} \\
	--output-reference Reference.matched.cnn \\
	--output-dir CNVKit/matched/ \\
	--short-names \\
	--diagram \\
	--scatter \\
	--annotate ${reference_flat} \\
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
		mv \${file} \${file_new}
done

	"""

	stub:
	"""#!/usr/bin/env bash
mkdir -p CNVKit/matched
#cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.matched.cns CNVKit/matched/
touch CNVKit/matched/${meta.sampleName}.matched.cns
	"""

}

process cnv_kit_single {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/CNVKit", mode: "copy", pattern: "CNVKit/single/${meta.sampleName}.${type}.cns", saveAs: { it.replaceFirst ("^CNVKit/single/","") }

	input:
		val (genome_build)
		val (reference)
		val (reference_flat)
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
mkdir -p CNVKit/matched
cnvkit.py batch \\
	${bam} \\
	--normal \\
	--fasta ${reference} \\
	--output-reference Reference.${type}.cnn \\
	--output-dir CNVKit/single/ \\
	--short-names \\
	--diagram \\
	--scatter \\
	--annotate ${reference_flat} \\
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
		mv \${file} \${file_new}
done

	"""

	stub:
	"""#!/usr/bin/env bash
mkdir -p CNVKit/single
#cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.${type}.cns CNVKit/single/
touch CNVKit/single/${meta.sampleName}.${type}.cns
	"""
}


process cnv_kit_segment {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/CNVKit", mode: "copy", pattern: "*.cns"

	input:
		val (genome_build)
		val (coverage_source)
		tuple val (meta), val (type), val(resolution), path (coverage_cnr)

	output:
		tuple val (meta), val (type), val ("${coverage_source}-cnv-kit"), val (resolution), path (coverage_cnr), path ("${meta.sampleName}.${type}.${coverage_source}.cns"), emit: result
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.${coverage_source}.cns"), emit: cns

	script:
	"""#!/usr/bin/env bash

cnvkit.py segment -o ${meta.sampleName}.${type}.${coverage_source}.cns ${coverage_cnr}
	"""

	stub:
	"""#!/usr/bin/env bash
#cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/CNVKit/${meta.sampleName}.${type}.${coverage_source}.cns .
touch ${meta.sampleName}.${type}.${coverage_source}.cns
	"""
}


