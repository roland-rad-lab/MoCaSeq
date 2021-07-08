
process loh_matched {
	tag "${meta.sampleName}"

	input:
		tuple path (interval_bed), path (interval_bed_index)
		tuple val (meta), path (normal_vcf), path (normal_vcf_index), path (tumor_vcf), path (tumor_vcf_index)

	output:
		tuple val (meta), path ("${meta.sampleName}.VariantsForLOHGermline.tsv.gz"), emit: result

	script:
	"""#!/usr/bin/env bash
source ${params.script_base}/file_handling.sh
temp_file_b=\$(moc_mktemp_file . bed)
trap "rm \${temp_file_b}" EXIT

extract_if_zip ${interval_bed} interval_bed_extracted \${temp_file_b}

bcftools merge \\
	--regions-file \${interval_bed_extracted} \\
	--merge all \\
	--output-type u \\
	${normal_vcf} \\
	${tumor_vcf} \\
	| bcftools view \\
	-m2 -M2 \\
	--types snps,indels \\
	--output-type u \\
	| bcftools filter \\
	--exclude 'INFO/MMQ[1]<60 & FORMAT/DP<10' \\
	--output-type u \\
	| bcftools query \\
	--format '[%SAMPLE\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%AF\\t%AD{0}\\t%AD{1}\\t%INFO/MMQ{1}\\t%INFO/MBQ{1}\\n]' \\
	| gzip > ${meta.sampleName}.VariantsForLOHGermline.tsv.gz
	"""
}

process loh_matched_assign_alleles {
	tag "${meta.sampleName}"

	input:
		tuple val (meta), path (loh_variants_tsv)

	script:
	"""#!/usr/bin/env python3
import gzip

with open ("${meta.sampleName}.VariantsForLOH.tsv", "w") as output_file:
	with gzip.open ("${loh_variants_tsv}", "rt") as input_file:
		for line in input_file:
			ldata = line.rstrip ().split ("\\t")
			print (ldata) 

	"""

}



