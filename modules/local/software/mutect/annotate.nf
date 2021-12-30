
params.gatk = [:]
params.mutect = [:]
params.snpeff = [:]
params.snpsift = [:]

process mutect_extract_single {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Mutect2", mode: "copy"

	input:
		val (genome_build)
		val (sift_fields)
		tuple val (meta), val (type), path (vcf), path (vcf_index)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.txt"), emit: result

	script:
	"""#!/usr/bin/env bash

zcat ${vcf} | ${params.snpeff.one_per_line} | bgzip -c > ${meta.sampleName}.${type}.one.vcf.gz

java -jar ${params.snpsift.jar} extractFields \\
	${meta.sampleName}.${type}.one.vcf.gz \\
	CHROM POS REF ALT "GEN[${type}].AF" "GEN[${type}].AD[0]" "GEN[${type}].AD[1]" \\
	${sift_fields} > ${meta.sampleName}.${type}.Mutect2.txt

	"""

	stub:
	"""#!/usr/bin/env bash

if [[ "${params.stub_json_map?.mutect_extract_single}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.Mutect2.txt .
fi

touch ${meta.sampleName}.${type}.Mutect2.txt

	"""

}

process mutect_extract_matched {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Mutect2", mode: "copy"

	input:
		val (genome_build)
		val (sift_fields)
		tuple val (meta), val(type), path (vcf), path (vcf_index)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.txt"), emit: result

	script:
	"""#!/usr/bin/env bash

zcat ${vcf} | ${params.snpeff.one_per_line} | bgzip -c > ${meta.sampleName}.${type}.one.vcf.gz

java -jar ${params.snpsift.jar} extractFields \\
	${meta.sampleName}.${type}.one.vcf.gz \\
	CHROM POS REF ALT "GEN[Tumor].AF" \\
	"GEN[Tumor].AD[0]" "GEN[Tumor].AD[1]" \\
	"GEN[Normal].AD[0]" "GEN[Normal].AD[1]" \\
	${sift_fields} > ${meta.sampleName}.${type}.Mutect2.txt

	"""

	stub:
	"""#!/usr/bin/env bash

if [[ "${params.stub_json_map?.mutect_extract_matched}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.Mutect2.txt .
fi

touch ${meta.sampleName}.${type}.Mutect2.txt

	"""

}

process mutect_filter_result_impact {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Mutect2", mode: "copy"

	input:
		val (genome_build)
		val (cgc_tsv)
		tuple val (meta), val (type), path (result_tsv)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.txt")
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.txt")
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt")

	script:
	"""#!/usr/bin/env Rscript

library (dplyr)

data <- read.table (file="${result_tsv}",sep="\\t",header=T,stringsAsFactors=F)
head (data)

data_cgc <- read.table (file="${cgc_tsv}",sep="\\t",header=T,stringsAsFactors=F)
head (data_cgc)

data_rare_impact <- data %>%
	dplyr::filter (ANN....IMPACT %in% c("HIGH", "MODERATE")) %>%
	data.frame

write.table (data_rare_impact %>% dplyr::mutate (dplyr::across (dplyr::everything (),~ ifelse (is.na (.x),"",.x))) %>% data.frame,file="${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.txt",sep="\\t",row.names=F,quote=F)

data_rare_impact_cgc <- data_rare_impact %>%
	dplyr::left_join (data_cgc %>% dplyr::select (Gene.Symbol) %>% data.frame,by=c("ANN....GENE"="Gene.Symbol")) %>%
	data.frame

write.table (data_rare_impact_cgc %>% dplyr::mutate (dplyr::across (dplyr::everything (),~ ifelse (is.na (.x),"",.x))) %>% data.frame,file="${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt",sep="\\t",row.names=F,quote=F)

	"""
}

process mutect_filter_result_impact_rare {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Mutect2", mode: "copy"

	input:
		val (genome_build)
		val (cgc_tsv)
		val (tru_sight)
		tuple val (meta), val (type), path (result_tsv)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.txt")
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.txt")
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt")
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.CGC.TruSight.txt")

	script:
	"""#!/usr/bin/env Rscript

library (dplyr)

data <- read.table (file="${result_tsv}",sep="\\t",header=T,stringsAsFactors=F)
head (data)

data_cgc <- read.table (file="${cgc_tsv}",sep="\\t",header=T,stringsAsFactors=F)
head (data_cgc)

data_tru_sight <- read.table (file="${tru_sight}",sep="\\t",header=F,stringsAsFactors=F)
head (data_tru_sight)

data_rare <- data %>%
	dplyr::mutate (dplyr::across (AF, ~ dplyr::if_else (is.na (.x),0,.x))) %>%
	dplyr::mutate (dplyr::across (c(AC,AN), ~ dplyr::if_else (is.na (.x),0L,.x))) %>%
	dplyr::filter (G5=="false",AF < 0.1 & AN < 100 | AF <0.01 & AN >= 100) %>%
	data.frame

write.table (data_rare %>% dplyr::mutate (dplyr::across (dplyr::everything (),~ ifelse (is.na (.x),"",.x))) %>% data.frame,file="${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.txt",sep="\\t",row.names=F,quote=F)

data_rare_impact <- data_rare %>%
	dplyr::filter (ANN....IMPACT %in% c("HIGH", "MODERATE")) %>%
	data.frame

write.table (data_rare_impact %>% dplyr::mutate (dplyr::across (dplyr::everything (),~ ifelse (is.na (.x),"",.x))) %>% data.frame,file="${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.txt",sep="\\t",row.names=F,quote=F)

data_rare_impact_cgc <- data_rare_impact %>%
	dplyr::left_join (data_cgc %>% dplyr::select (Gene.Symbol) %>% data.frame,by=c("ANN....GENE"="Gene.Symbol")) %>%
	data.frame

write.table (data_rare_impact_cgc %>% dplyr::mutate (dplyr::across (dplyr::everything (),~ ifelse (is.na (.x),"",.x))) %>% data.frame,file="${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt",sep="\\t",row.names=F,quote=F)

data_rare_impact_cgc_tru_sight <- data_rare_impact_cgc %>%
	dplyr::left_join (data_tru_sight %>% dplyr::select (V1) %>% data.frame,by=c("ANN....GENE"="V1")) %>%
	data.frame

write.table (data_rare_impact_cgc_tru_sight %>% dplyr::mutate (dplyr::across (dplyr::everything (),~ ifelse (is.na (.x),"",.x))) %>% data.frame,file="${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.CGC.TruSight.txt",sep="\\t",row.names=F,quote=F)

	"""

}

process mutect_filter {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Mutect2", mode: "copy"

	input:
		val (genome_build)
		val (reference)
		tuple val (meta), val (type), path (vcf), path (vcf_index), path (stats), path (model)

	output:
		tuple val (meta), val(type), path ("${meta.sampleName}.${type}.m2.filtered.vcf.gz"), path ("${meta.sampleName}.${type}.m2.filtered.vcf.gz.tbi"), emit: result
		path ("${meta.sampleName}.${type}.m2.filtered.summary.txt")

	script:
	"""#!/usr/bin/env bash

if [[ "${params.mutect.artefact}" == "no" ]]; then
	java -jar ${params.gatk.jar} FilterMutectCalls \\
	--variant ${vcf} \\
	--output ${meta.sampleName}.${type}.m2.filtered.vcf.gz \\
	--reference ${reference}
elif [[ "${params.mutect.artefact}" == "yes" ]]; then
	java -jar ${params.gatk.jar} FilterMutectCalls \\
	--variant ${vcf} \\
	--output ${meta.sampleName}.${type}.m2.filtered.vcf.gz \\
	--reference ${reference} \\
	--ob-priors ${model}
else
	echo "Please supply a valid value for mutect.artefact '${params.mutect.artefact}' must be one of yes,no"
	exit 1
fi


# output filtering statistics
zcat ${meta.sampleName}.${type}.m2.filtered.vcf.gz | grep "^[^#;]" | cut -f 7 | sort | uniq -c | sort -nr > ${meta.sampleName}.${type}.m2.filtered.summary.txt

	"""

	stub:
	"""#!/usr/bin/env bash

if [[ "${params.stub_json_map?.mutect_filter}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.m2.filtered.vcf.gz .
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.m2.filtered.vcf.gz.tbi .
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.m2.filtered.summary.txt .
fi

touch ${meta.sampleName}.${type}.m2.filtered.vcf.gz
touch ${meta.sampleName}.${type}.m2.filtered.vcf.gz.tbi
touch ${meta.sampleName}.${type}.m2.filtered.summary.txt

	"""

}

process mutect_post_process_single
{
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Mutect2", mode: "copy"

	input:
		val (genome_build)
		val (all_snp_file)
		val (common_snp_file)
		tuple val (meta), val (type), path (vcf), path (vcf_index)

	output:
		tuple val (meta), val(type), path ("${meta.sampleName}.${type}.Mutect2.vcf.gz"), path ("${meta.sampleName}.${type}.Mutect2.vcf.gz.tbi"), emit: result
		path ("${meta.sampleName}.${type}.Mutect2.Positions.txt")

	script:
	"""#!/usr/bin/env bash

java -jar ${params.snpsift.jar} extractFields \\
${vcf} \\
CHROM POS REF ALT "GEN["$type"].AF" "GEN["$type"].AD[0]" \\
"GEN["$type"].AD[1]" MMQ[1] MBQ[1] \\
> ${meta.sampleName}.${type}.Mutect2.Positions.txt

java -jar ${params.gatk.jar} SelectVariants --max-indel-size 10 \\
-V ${vcf} \\
-output ${meta.sampleName}.${type}.m2.filtered.selected.vcf.gz

if [[ "${params.mutect.filter}" = "soft" ]]; then
	zcat ${meta.sampleName}.${type}.m2.filtered.selected.vcf.gz \\
	| java -jar ${params.snpsift.jar} filter \\
	"( ( FILTER = 'PASS') & (GEN[${type}].AF >= 0.05) & \\
	(GEN[${type}].AD[1] >= 2) & (GEN[${type}].AD[0] + GEN[${type}].AD[1] >= 5) )" \\
	| bgzip -c > ${meta.sampleName}.${type}.m2.postprocessed.vcf.gz

	tabix -p vcf ${meta.sampleName}.${type}.m2.postprocessed.vcf.gz

	bcftools isec -C -c none -O z -w 1 \\
	-o ${meta.sampleName}.${type}.m2.postprocessed.snp_removed.vcf.gz \\
	${meta.sampleName}.${type}.m2.postprocessed.vcf.gz \\
	${common_snp_file}

elif [[ "${params.mutect.filter}" = "hard" ]]; then
	zcat ${meta.sampleName}.${type}.m2.filtered.selected.vcf.gz \\
	| java -jar ${params.snpsift.jar} filter \\
	"( ( FILTER = 'PASS') & (GEN[${type}].AF >= 0.1) & \\
	(GEN[${type}].AD[1] >= 2) & (GEN[${type}].AD[0] + GEN[${type}].AD[1] >= 10) )" \\
	| bgip -c > ${meta.sampleName}.${type}.m2.postprocessed.vcf.gz

	tabix -p vcf ${meta.sampleName}.${type}.m2.postprocessed.vcf.gz

	bcftools isec -C -c none -O z -w 1 \\
	-o ${meta.sampleName}.${type}.m2.postprocessed.snp_removed.vcf.gz \\
	${meta.sampleName}.${type}.m2.postprocessed.vcf.gz \\
	${all_snp_file}

elif [[ "${params.mutect.filter}" = "none" ]]; then
	zcat ${meta.sampleName}.${type}.m2.filtered.selected.vcf.gz \\
	| java -jar ${params.snpsift.jar} filter \\
	"( ( FILTER = 'PASS' ) )" \\
	| bgzip -c > ${meta.sampleName}.${type}.m2.postprocessed.snp_removed.vcf.gz

else
	echo "Please supply a valid value for mutect.filter '${params.mutect.filter}' must be one of soft,hard,none"
	exit 1
fi

bcftools norm -m -any -O z \\
	-o ${meta.sampleName}.${type}.Mutect2.vcf.gz \\
	${meta.sampleName}.${type}.m2.postprocessed.snp_removed.vcf.gz

tabix -p vcf ${meta.sampleName}.${type}.Mutect2.vcf.gz

	"""

	stub:
	"""#!/usr/bin/env bash

if [[ "${params.stub_json_map?.mutect_post_process_single}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.Mutect2.vcf.gz .
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.Mutect2.vcf.gz.tbi .
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.Mutect2.Positions.txt .
fi

touch ${meta.sampleName}.${type}.Mutect2.vcf.gz
touch ${meta.sampleName}.${type}.Mutect2.vcf.gz.tbi
touch ${meta.sampleName}.${type}.Mutect2.Positions.txt

	"""

}

process mutect_post_process_matched {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Mutect2", mode: "copy"

	input:
		val (genome_build)
		val (all_snp_file)
		val (common_snp_file)
		tuple val (meta), val (type), path (vcf), path (vcf_index)

	output:
		tuple val (meta), val(type), path ("${meta.sampleName}.${type}.Mutect2.vcf.gz"), path ("${meta.sampleName}.${type}.Mutect2.vcf.gz.tbi"), emit: result

	script:
	"""#!/usr/bin/env bash

java -jar ${params.gatk.jar} SelectVariants \\
	--max-indel-size 10 \\
	-V ${vcf} \\
	-output ${meta.sampleName}.${type}.m2.filtered.selected.vcf.gz

if [[ "${params.mutect.filter}" == "soft" ]]; then
	zcat ${meta.sampleName}.${type}.m2.filtered.selected.vcf.gz \\
	| java -jar ${params.snpsift.jar} filter \\
	"( ( FILTER = 'PASS') & (GEN[Tumor].AF >= 0.05) & \\
	( ( GEN[Tumor].AD[0] + GEN[Tumor].AD[1]) >= 5 ) & \\
	( ( GEN[Normal].AD[0] + GEN[Normal].AD[1]) >= 5 ) & \\
	(GEN[Tumor].AD[1] >= 2) & (GEN[Normal].AD[1] <= 1) )" \\
	| bgzip -c > ${meta.sampleName}.${type}.m2.postprocessed.vcf.gz

	tabix -p vcf ${meta.sampleName}.${type}.m2.postprocessed.vcf.gz

	bcftools isec -C -c none -O z -w 1 \\
	-o ${meta.sampleName}.${type}.m2.postprocessed.snp_removed.vcf.gz \\
	${meta.sampleName}.${type}.m2.postprocessed.vcf.gz \\
	${common_snp_file}

elif [[ "${params.mutect.filter}" == "hard" ]]; then
	zcat ${meta.sampleName}.${type}.m2.filtered.selected.vcf.gz \\
	| java -jar ${params.snpsift.jar} filter \\
	"( ( FILTER = 'PASS') & (GEN[Tumor].AF >= 0.1) & \\
	( ( GEN[Tumor].AD[0] + GEN[Tumor].AD[1]) >= 10 ) & \\
	( ( GEN[Normal].AD[0] + GEN[Normal].AD[1]) >= 10 ) & \\
	(GEN[Tumor].AD[1] >= 3) & (GEN[Normal].AD[1] = 0) )" \\
	| bgzip -c > ${meta.sampleName}.${type}.m2.postprocessed.vcf.gz

	tabix -p vcf ${meta.sampleName}.${type}.m2.postprocessed.vcf.gz

	bcftools isec -C -c none -O z -w 1 \\
	-o ${meta.sampleName}.${type}.m2.postprocessed.snp_removed.vcf.gz \\
	${meta.sampleName}.${type}.m2.postprocessed.vcf.gz \\
	${all_snp_file}

elif [[ "${params.mutect.filter}" == "none" ]]; then
	zcat ${meta.sampleName}.${type}.m2.filtered.selected.vcf.gz \\
	| java -jar ${params.snpsift.jar} filter \\
	"( ( FILTER = 'PASS' ) )" \\
	| bgzip -c > ${meta.sampleName}.${type}.m2.postprocessed.snp_removed.vcf.gz

	tabix -p vcf ${meta.sampleName}.${type}.m2.postprocessed.snp_removed.vcf.gz
else
	echo "Please supply a valid value for mutect.filter '${params.mutect.filter}' must be one of soft,hard,none"
	exit 1
fi

bcftools norm -m -any -O z \\
	-o ${meta.sampleName}.${type}.Mutect2.vcf.gz \\
	${meta.sampleName}.${type}.m2.postprocessed.snp_removed.vcf.gz

tabix -p vcf ${meta.sampleName}.${type}.Mutect2.vcf.gz

	"""

	stub:
	"""#!/usr/bin/env bash

if [[ "${params.stub_json_map?.mutect_post_process_single}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.Mutect2.vcf.gz .
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.Mutect2.vcf.gz.tbi .
fi

touch ${meta.sampleName}.${type}.Mutect2.vcf.gz
touch ${meta.sampleName}.${type}.Mutect2.vcf.gz.tbi

	"""

}

process mutect_sift {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Mutect2", mode: "copy"

	input:
		val (genome_build)
		val (dbnsfp)
		val (sift_sources)
		val (snpeff_version)
		tuple val (meta), val (type), path (vcf), path (vcf_index)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.annotated.vcf.gz"), path ("${meta.sampleName}.${type}.Mutect2.annotated.vcf.gz.tbi"), emit: result
		path ("${meta.sampleName}.${type}.Mutect2.annotated.vcf.gz.stats")

	script:
	switch ( meta["organism"] )
	{
		case "human":
			"""#!/usr/bin/env bash

ln -s \$(readlink ${vcf}) input.0.vcf.gz
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

ticker=\$((ticker_prev+1))
java -Xmx16g -jar ${params.snpsift.jar} \\
	DbNSFP -db ${dbnsfp} \\
	-v input.\${ticker_prev}.vcf.gz \\
	-f MetaLR_pred,MetaSVM_pred,SIFT_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,PROVEAN_pred \\
	| bgzip -c > input.\${ticker}.vcf.gz
ticker_prev=\${ticker}

java -Xmx16g -jar ${params.snpeff.jar} ${snpeff_version} -canon \\
	-csvStats ${meta.sampleName}.${type}.Mutect2.annotated.vcf.gz.stats \\
	input.\${ticker_prev}.vcf.gz \\
	| bgzip -c > ${meta.sampleName}.${type}.Mutect2.annotated.vcf.gz

tabix -p vcf ${meta.sampleName}.${type}.Mutect2.annotated.vcf.gz

			"""
			break
		case "mouse":
			"""#!/usr/bin/env bash
ln -s \$(readlink ${vcf}) ${meta.sampleName}.${type}.Mutect2.annotated.vcf.gz
ln -s \$(readlink ${vcf_index}) ${meta.sampleName}.${type}.Mutect2.annotated.vcf.gz.tbi
			"""
			break
		default:
			exit 1, "[MoCaSeq] Error: Unrecogised organism '${meta.organism}' must be one of [human,mouse] in mutect_sift"
			break
	}

	stub:
	"""#!/usr/bin/env bash

if [[ "${params.stub_json_map?.mutect_sift}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.Mutect2.annotated.vcf.gz .
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.Mutect2.annotated.vcf.gz.tbi .
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.Mutect2.annotated.vcf.gz.stats .
fi

touch ${meta.sampleName}.${type}.Mutect2.annotated.vcf.gz
touch ${meta.sampleName}.${type}.Mutect2.annotated.vcf.gz.tbi
touch ${meta.sampleName}.${type}.Mutect2.annotated.vcf.gz.stats

	"""
}


