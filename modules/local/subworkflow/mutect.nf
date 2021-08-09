#!/usr/bin/env nextflow

include {
	mutect_matched;
	mutect_single as mutect_single_normal;
	mutect_single as mutect_single_tumor;
	mutect_combine_vcf;
	mutect_combine_vcf as mutect_combine_vcf_single_normal;
	mutect_combine_vcf as mutect_combine_vcf_single_tumor;
} from "../software/mutect/main"

workflow MUTECT
{
	take:
		genome_build
		ch_fasta
		ch_interval
		interval_n
		ch_data
	main:
		ch_data_expanded = ch_data.map { it ->
			tuple ( it, it["normalBAM"], it["tumorBAM"] )
		}

		ch_data_expanded_single_normal = ch_data.map { it ->
			tuple ( it, "Normal", it["normalBAM"] )
		}

		ch_data_expanded_single_tumor = ch_data.map { it ->
			tuple (it, "Tumor", it["tumorBAM"] )
		}

		mutect_matched (genome_build, ch_fasta, ch_data_expanded, ch_interval)
		mutect_single_normal (genome_build, ch_fasta, ch_data_expanded_single_normal, ch_interval)
		mutect_single_tumor (genome_build, ch_fasta, ch_data_expanded_single_tumor, ch_interval)

		ch_vcf = mutect_matched.out.result.map { [[it[1], it[0]["sampleName"]].join ("__"), it] }
			.groupTuple (size: interval_n.value)
			.map { it[1] }
			.map { [ it[0][0], it[0][1], it.collect { jt -> jt[2] }, it.collect { jt -> jt[3] } ] }

		ch_vcf_single_normal = mutect_single_normal.out.result.map { [[it[1], it[0]["sampleName"]].join ("__"), it] }
			.groupTuple (size: interval_n.value)
			.map { it[1] }
			.map { [ it[0][0], it[0][1], it.collect { jt -> jt[2] }, it.collect { jt -> jt[3] } ] }

		ch_vcf_single_tumor = mutect_single_tumor.out.result.map { [[it[1], it[0]["sampleName"]].join ("__"), it] }
			.groupTuple (size: interval_n.value)
			.map { it[1] }
			.map { [ it[0][0], it[0][1], it.collect { jt -> jt[2] }, it.collect { jt -> jt[3] } ] }

		mutect_combine_vcf (genome_build, ch_vcf)
		mutect_combine_vcf_single_normal (genome_build, ch_vcf_single_normal)
		mutect_combine_vcf_single_tumor (genome_build, ch_vcf_single_tumor)
	emit:
		result = mutect_combine_vcf.out.result.mix (mutect_combine_vcf_single_normal.out.result,mutect_combine_vcf_single_tumor.out.result)
}

