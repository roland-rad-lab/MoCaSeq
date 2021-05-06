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
		println interval_n
		println interval_n.value

		ch_vcf = mutect_matched.out.vcf.map { it ->
				[["matched", it[0]["sampleName"]].join ("__"), it]
			}
			.view ()
			.groupTuple (size: 24)
			.map { it -> it[1] }

		mutect_combine_vcf (ch_vcf)
	emit:
		result = mutect_combine_vcf.out.result
}

