
params.strelka = [:]

process strelka_matched {
	tag "${meta.sampleName}"

	input:
		val (reference)
		tuple val (interval_bed), val (interval_bed_index)
		tuple val (meta), path (bam_normal), path (bai_normal), path (bam_tumor), path (bai_tumor), path (candidate_small_indels_vcf), path (candidate_small_indels_vcf_index)

	output:
		tuple val (meta), val("matched"), path("Strelka/results/variants/somatic.snvs.vcf.gz"), path ("Strelka/results/variants/somatic.snvs.vcf.gz.tbi"), path ("Strelka/results/variants/somatic.indels.vcf.gz"), path ("Strelka/results/variants/somatic.indels.vcf.gz.tbi"), emit: result

	script:
	"""#!/usr/bin/env bash
python2 ${params.strelka.dir}/bin/configureStrelkaSomaticWorkflow.py \\
	--normalBam ${bam_normal} \\
	--tumorBam ${bam_tumor} \\
	--ref ${reference} \\
	--runDir Strelka \\
	--indelCandidates ${candidate_small_indels_vcf} \\
	--callRegions ${interval_bed}

python2 Strelka/runWorkflow.py -m local -j ${params.strelka.threads}
	"""
}

process strelka_matched_post {
	tag "${meta.sampleName}"

	input:
		tuple val (meta), val (type), path (somatic_snv_vcf), path (somatic_snv_vcf_index), path (somatic_indel_vcf), path (somatic_indel_vcf_index)

	output:
		tuple val (meta), path ("${meta.sampleName}.str.snp.post-processed.vcf.gz"), path ("${meta.sampleName}.str.snp.post-processed.vcf.gz.tbi"), path ("${meta.sampleName}.str.indel.post-processed.vcf.gz"), path ("${meta.sampleName}.str.indel.post-processed.vcf.gz.tbi"), emit: result

	script:
	"""#!/usr/bin/env bash

zcat somatic.snvs.vcf \\
	| java -jar ${params.annotation.snpeff.dir}/SnpSift.jar filter \\
	" \\
	( ( ( FILTER = 'PASS' ) & ( ALT = 'T' ) & ( GEN[NORMAL].TU[0] <= 1 ) & ( GEN[TUMOR].TU[0] >= 2 ) & ( GEN[TUMOR].DP >= 5 ) & ( GEN[NORMAL].DP >= 5 ) & ( GEN[TUMOR].DP * 0.05 <= GEN[TUMOR].TU[0] ) ) \\
	| ( ( FILTER = 'PASS' ) & ( ALT = 'A' ) & ( GEN[NORMAL].AU[0] <= 1 ) & ( GEN[TUMOR].AU[0] >= 2 ) & ( GEN[TUMOR].DP >= 5 ) & ( GEN[NORMAL].DP >= 5 ) & ( GEN[TUMOR].DP * 0.05 <= GEN[TUMOR].AU[0] ) ) \\
	| ( ( FILTER = 'PASS' ) & ( ALT = 'C' ) & ( GEN[NORMAL].CU[0]<= 1 ) & ( GEN[TUMOR].CU[0] >= 2 ) & ( GEN[TUMOR].DP >= 5 ) & ( GEN[NORMAL].DP >= 5 ) & ( GEN[TUMOR].DP * 0.05 <= GEN[TUMOR].CU[0] ) ) \\
	| ( ( FILTER = 'PASS' ) & ( ALT = 'G' ) & ( GEN[NORMAL].GU[0] <= 1 ) & ( GEN[TUMOR].GU[0] >= 2 ) & ( GEN[TUMOR].DP >= 5 ) & ( GEN[NORMAL].DP >= 5 ) & ( GEN[TUMOR].DP * 0.05 <= GEN[TUMOR].GU[0] ) ) ) \\
	" \\
	| bgzip -c > ${meta.sampleName}.str.snp.post-processed.vcf.gz
tabix -p vcf ${meta.sampleName}.str.snp.post-processed.vcf.gz

zcat ${somatic_indel_vcf} \\
	| java -jar ${params.annotation.snpeff.dir}/SnpSift.jar filter \\
	" \\
	( ( FILTER = 'PASS' ) & ( GEN[TUMOR].DP >= 5 ) & ( GEN[NORMAL].DP >= 5 ) & ( GEN[NORMAL].TIR[0] <= 1) & \\
	( GEN[TUMOR].TIR[0] >= 2 ) & ( GEN[TUMOR].DP * 0.05 <= GEN[TUMOR].TIR[0] ) )
	" \\
	> ${meta.sampleName}.str.indel._post-processed.vcf

java -jar ${params.gatk.jar} SelectVariants \\
	--max-indel-size 10 \\
	--variant ${meta.sampleName}.str.indel._post-processed.vcf \\
	--output ${meta.sampleName}.str.indel.post-processed.vcf.gz

	"""
}



