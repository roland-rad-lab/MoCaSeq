#!/usr/bin/env nextflow

include { mutect_matched; mutect_combine_vcf } from "../software/mutect/main"

workflow MUTECT
{
	take:
		ch_fasta
		ch_interval
		interval_n
		ch_data
	main:
		ch_data_expanded = ch_data.map { it ->
			tuple ( it, it["normalBAM"], it["tumorBAM"] )
		}

		mutect_matched (ch_fasta, ch_data_expanded, ch_interval)

		ch_vcf = mutect_matched.out.vcf.map { [["matched", it[0]["sampleName"]].join ("__"), it] }
			.groupTuple (size: interval_n.value)
			.map { it[1] }
			.map { [ it[0][0], it.collect { jt -> jt[1] } ] }

		mutect_combine_vcf (ch_vcf)
	emit:
		result = mutect_combine_vcf.out.result
}

