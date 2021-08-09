#!/usr/bin/env nextflow

include { manta_matched; manta_matched_post } from "../software/manta/main"

workflow MANTA
{
	take:
		genome_build
		ch_fasta
		ch_interval_bed
		ch_data

	main:
		ch_data_expanded = ch_data.map { it ->
			tuple ( it, it["normalBAM"], it["normalBAI"], it["tumorBAM"], it["tumorBAI"] )
		}

		manta_matched (ch_fasta, ch_interval_bed, ch_data_expanded)
		manta_matched_post (manta_matched.out.sv)

	emit:
		result = manta_matched_post.out.result
		basic = manta_matched_post.out.filtered_vcf
		indel = manta_matched.out.indel
}

