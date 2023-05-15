#!/usr/bin/env nextflow

include {
	mutect_matched;
	mutect_single as mutect_single_normal;
	mutect_single as mutect_single_tumor;
	mutect_combine_vcf;
	mutect_combine_vcf as mutect_combine_vcf_single_normal;
	mutect_combine_vcf as mutect_combine_vcf_single_tumor;
} from "../software/mutect/main"

include {
	mutect_filter;
	mutect_extract_single;
	mutect_extract_matched;
	mutect_post_process_single;
	mutect_post_process_matched;
	mutect_sift
} from "../software/mutect/annotate"

include {
	mutect_filter_result_impact;
	mutect_filter_result_impact_rare;
	mutect_signatures_matched
} from "../software/mutect/result"

workflow MUTECT
{
	take:
		genome_build
		ch_fasta
		ch_interval
		interval_n
		ch_data
	main:
		if (params.debug) { println "[MoCaSeq] debug: entered MUTECT subworkflow" }
		ch_data_expanded = ch_data.filter { it["type"] == "Tumor" }.map { it ->
			tuple ( it, it["normalBAM"], it["normalBAI"], it["tumorBAM"], it["tumorBAI"] )
		}

		ch_data_single_branched = ch_data.branch {
			normal: it["type"] == "Normal"
			tumor: it["type"] = "Tumor"
			other: true
		}

		ch_data_single_branched.other.view { "[MoCaSeq] error: Unknown (type) for input:\n${it}\nExpected: [Normal,Tumor]." }

		mutect_matched (genome_build, ch_fasta, ch_data_expanded, ch_interval)
		if (params.debug) { 
			println "[MoCaSeq] debug: pre mutect_single_normal process"
			ch_data_single_branched.normal.view()
			ch_data_expanded_normal.view()
		}
		mutect_single_normal (genome_build, ch_fasta, ch_data_single_branched.normal.map { it -> tuple ( it, "Normal", it["normalBAM"], it["normalBAI"] ) }, ch_interval)
		if (params.debug) { 
			println "[MoCaSeq] debug: pre mutect_single_tumor process"
			ch_data_single_branched.tumor.view()
			ch_data_expanded_normal.view()
		}
		mutect_single_tumor (genome_build, ch_fasta, ch_data_single_branched.tumor.map { it -> tuple (it, "Tumor", it["tumorBAM"], it["tumorBAI"] ) }, ch_interval)

		ch_vcf = mutect_matched.out.result.map { [[it[1], it[0]["sampleName"]].join ("__"), it] }
			.groupTuple (size: interval_n.value)
			.map { it[1] }
			.map { [ it[0][0], it[0][1], it.collect { jt -> jt[2] }, it.collect { jt -> jt[3] }, it.collect { jt -> jt[4] } ] }

		ch_vcf_single_normal = mutect_single_normal.out.result.map { [[it[1], it[0]["sampleName"]].join ("__"), it] }
			.groupTuple (size: interval_n.value)
			.map { it[1] }
			.map { [ it[0][0], it[0][1], it.collect { jt -> jt[2] }, it.collect { jt -> jt[3] }, it.collect { jt -> jt[4] } ] }

		ch_vcf_single_tumor = mutect_single_tumor.out.result.map { [[it[1], it[0]["sampleName"]].join ("__"), it] }
			.groupTuple (size: interval_n.value)
			.map { it[1] }
			.map { [ it[0][0], it[0][1], it.collect { jt -> jt[2] }, it.collect { jt -> jt[3] }, it.collect { jt -> jt[4] } ] }

		mutect_combine_vcf (genome_build, ch_vcf)
		mutect_combine_vcf_single_normal (genome_build, ch_vcf_single_normal)
		mutect_combine_vcf_single_tumor (genome_build, ch_vcf_single_tumor)
	emit:
		result = mutect_combine_vcf.out.result.mix (mutect_combine_vcf_single_normal.out.result,mutect_combine_vcf_single_tumor.out.result)
		full = mutect_combine_vcf.out.full.mix (mutect_combine_vcf_single_normal.out.full,mutect_combine_vcf_single_tumor.out.full)
}

workflow MUTECT_ANNOTATE
{
	take:
		genome_build
		ch_fasta
		ch_data
		ch_snpeff_version
		ch_all_vcf
		ch_common_vcf
		ch_dbnsfp
		ch_sift_sources
		ch_sift_fields

	main:
		mutect_filter (genome_build, ch_fasta, ch_data)

		ch_filter_branched = mutect_filter.out.result.branch {
			single: it[1] == "Tumor" || it[1] == "Normal"
			matched: it[1] == "matched"
			other: true
		}

		ch_filter_branched.other.view { "[MoCaSeq] error: Unknown (type) for input:\n${it}\nExpected: [Tumor,Normal,matched]." }

		mutect_post_process_single (genome_build, ch_all_vcf, ch_common_vcf, ch_filter_branched.single)
		mutect_post_process_matched (genome_build, ch_all_vcf, ch_common_vcf, ch_filter_branched.matched)

		mutect_sift (genome_build, ch_dbnsfp, ch_sift_sources, ch_snpeff_version, mutect_post_process_single.out.result.mix (mutect_post_process_matched.out.result))

		ch_sift_branched = mutect_sift.out.result.branch {
			single: it[1] == "Tumor" || it[1] == "Normal"
			matched: it[1] == "matched"
			other: true
		}

		ch_sift_branched.other.view { "[MoCaSeq] error: Unknown (type) for input:\n${it}\nExpected: [Tumor,Normal,matched]." }

		mutect_extract_single (genome_build, ch_sift_fields, ch_sift_branched.single)
		mutect_extract_matched (genome_build, ch_sift_fields, ch_sift_branched.matched)

	emit:
		result = mutect_extract_single.out.result.mix (mutect_extract_matched.out.result)
		post_process = mutect_post_process_single.out.result.mix (mutect_post_process_matched.out.result)
}

workflow MUTECT_RESULT
{
	take:
		genome_build
		ch_data_post_process
		ch_data
		ch_cgc
	main:
		mutect_filter_result_impact (genome_build, ch_cgc, ch_data)
		mutect_signatures_matched (genome_build, ch_data_post_process.filter { it[1] == "matched" })
}

workflow MUTECT_RESULT_RARE
{

	take:
		genome_build
		ch_data_post_process
		ch_data
		ch_cgc
		ch_tru_sight

	main:
		mutect_filter_result_impact_rare (genome_build, ch_cgc, ch_tru_sight, ch_data)
		mutect_signatures_matched (genome_build, ch_data_post_process.filter { it[1] == "matched" })
}


