
params.gatk = [:]
params.mutect = [:]
params.snpeff = [:]
params.snpsift = [:]

process mutect_extract {
	tag "${meta.sampleName}"

	input:
		val (genome_build)
		tuple val (meta), path (vcf)
		val (sift_fields)

	script:
	"""#!/usr/bin/env bash

	"""
}

process mutect_filter {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Mutect2", mode: "copy"

	input:
		val (genome_build)
		val (reference)
		val (all_snp_file)
		val (common_snp_file)
		tuple val (meta), val (type), path (vcf), path (vcf_index), path (model)

	output:
		tuple val (meta), path ("${meta.sampleName}.Mutect2.vcf.gz"), path ("${meta.sampleName}.Mutect2.vcf.gz.tbi"), emit: result
		path ("${meta.sampleName}.Mutect.filtered.summary.txt")

	script:
	"""#!/usr/bin/env bash

if [[ "${params.mutect.artefact}" == "no" ]]; then
	# just copy without filtering
	java -jar ${params.gatk.jar} FilterMutectCalls \\
	--variant ${vcf} \\
	--output ${meta.sampleName}.Mutect.filtered.vcf.gz \\
	--reference ${reference}

	#cp ${vcf} ${meta.sampleName}.Mutect.filtered.vcf.gz
	#cp ${vcf_index} ${meta.sampleName}.Mutect.filtered.vcf.gz.tbi
elif [[ "${params.mutect.artefact}" == "yes" ]]; then
	# We expect the ob-file to already be available
	# filter artifacts
	java -jar ${params.gatk.jar} FilterMutectCalls \\
	--variant ${vcf} \\
	--output ${meta.sampleName}.Mutect.filtered.vcf.gz \\
	--reference ${reference} \\
	--ob-priors ${model}
fi

# output filtering statistics
zcat ${meta.sampleName}.Mutect.filtered.vcf.gz | grep "^[^#;]" | cut -f 7 | sort | uniq -c | sort -nr > ${meta.sampleName}.Mutect.filtered.summary.txt

java -jar ${params.gatk.jar} SelectVariants \\
	--max-indel-size 10 \\
	-V ${meta.sampleName}.Mutect.filtered.vcf.gz \\
	-output ${meta.sampleName}.Mutect.filtered.selected.vcf.gz

if [[ "${params.mutect.filter}" == "soft" ]]; then
	zcat ${meta.sampleName}.Mutect.filtered.selected.vcf.gz \\
	| java -jar ${params.snpsift.jar} filter \\
	"( ( FILTER = 'PASS') & (GEN[Tumor].AF >= 0.05) & \\
	( ( GEN[Tumor].AD[0] + GEN[Tumor].AD[1]) >= 5 ) & \\
	( ( GEN[Normal].AD[0] + GEN[Normal].AD[1]) >= 5 ) & \\
	(GEN[Tumor].AD[1] >= 2) & (GEN[Normal].AD[1] <= 1) )" \\
	| bgzip -c > ${meta.sampleName}.m2.postprocessed.vcf.gz

	tabix -p vcf ${meta.sampleName}.m2.postprocessed.vcf.gz

	bcftools isec -C -c none -O z -w 1 \\
	-o ${meta.sampleName}.m2.postprocessed.snp_removed.vcf.gz \\
	${meta.sampleName}.m2.postprocessed.vcf.gz \\
	${common_snp_file}

elif [[ "${params.mutect.filter}" == "hard" ]]; then
	zcat ${meta.sampleName}.Mutect.filtered.selected.vcf.gz \\
	| java -jar ${params.snpsift.jar} filter \\
	"( ( FILTER = 'PASS') & (GEN[Tumor].AF >= 0.1) & \\
	( ( GEN[Tumor].AD[0] + GEN[Tumor].AD[1]) >= 10 ) & \\
	( ( GEN[Normal].AD[0] + GEN[Normal].AD[1]) >= 10 ) & \\
	(GEN[Tumor].AD[1] >= 3) & (GEN[Normal].AD[1] = 0) )" \\
	| bgzip -c > ${meta.sampleName}.m2.postprocessed.vcf.gz

	tabix -p vcf ${meta.sampleName}.m2.postprocessed.vcf.gz

	bcftools isec -C -c none -O z -w 1 \\
	-o ${meta.sampleName}.m2.postprocessed.snp_removed.vcf.gz \\
	${meta.sampleName}.m2.postprocessed.vcf.gz \\
	${all_snp_file}

elif [[ "${params.mutect.filter}" == "none" ]]; then
	cat ${meta.sampleName}.Mutect.filtered.selected.vcf \\
	| java -jar ${params.snpsift.jar} filter \\
	"( ( FILTER = 'PASS' ) )" \\
	| bgzip -c > ${meta.sampleName}.m2.postprocessed.snp_removed.vcf.gz

	tabix -p vcf ${meta.sampleName}.m2.postprocessed.snp_removed.vcf.gz
else
	echo "Please supply a valid value for mutect.filter '${params.mutect.filter}' must be one of soft,hard,none"
	exit 1
fi

bcftools norm -m -any -O z \\
	-o zcat ${meta.sampleName}.Mutect2.vcf.gz \\
	${meta.sampleName}.m2.postprocessed.snp_removed.vcf.gz

tabix -p vcf ${meta.sampleName}.Mutect2.vcf.gz

	"""

}

process mutect_sift {
	tag "${meta.sampleName}"

	input:
		val (genome_build)
		tuple val (meta), path (vcf), path (vcf_index)
		val (sift_sources)

	output:
		tuple val (meta), path ("${meta.sampleName}.Mutect.filtered.annotated.vcf.gz"), path ("${meta.sampleName}.Mutect.filtered.annotated.vcf.gz.tbi"), emit: result

	script:
	"""#!/usr/bin/env bash

mv ${vcf} input.0.vcf.gz
ticker=1
ticker_prev=0
for sift_source in ${sift_sources};
do
	ticker=\$((ticker_prev+1))
	java -Xmx16g -jar ${params.snpsift.jar} annotate \\
	\${sift_source} input.\${ticker_prev}.vcf.gz \\
	| bgzip -c > input.\${ticker}.vcf.gz
	ticker_prev=\${ticker}
done

mv input.\${ticker_prev}.vcf.gz ${meta.sampleName}.Mutect.filtered.annotated.vcf.gz
tabix -p vcf ${meta.sampleName}.Mutect.filtered.annotated.vcf.gz

	"""
}


