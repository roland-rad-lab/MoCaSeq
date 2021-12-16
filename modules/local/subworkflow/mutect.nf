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
	mutect_extract;
	mutect_sift;
} from "../software/mutect/annotate"

workflow MUTECT
{
	take:
		genome_build
		ch_fasta
		ch_interval
		interval_n
		ch_data
	main:
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
		mutect_single_normal (genome_build, ch_fasta, ch_data_single_branched.normal.map { it -> tuple ( it, "Normal", it["normalBAM"], it["normalBAI"] ) }, ch_interval)
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
		ch_sift_sources
		ch_sift_fields

	main:
		mutect_filter (genome_build, ch_fasta, ch_snpeff_version, ch_all_vcf, ch_common_vcf, ch_data.filter { it[1] == "matched" })

		ch_filter_branched = mutect_filter.out.result.branch {
			human: it[0]["organism"] == "human"
			mouse: it[0]["organism"] == "mouse"
			other: true
		}

		ch_filter_branched.other.view { "[MoCaSeq] error: Unknown (type) for input:\n${it}\nExpected: [human,mouse]." }

		mutect_sift (genome_build, ch_filter_branched.human, ch_sift_sources)
		mutect_extract (genome_build, ch_filter_branched.mouse.mix (mutect_sift.out.result), ch_sift_fields)

}

