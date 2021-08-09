#!/usr/bin/env nextflow

include { delly_matched_call; delly_matched_filter } from "../software/delly/main"

workflow DELLY
{
	take:
		genome_build
		ch_fasta
		ch_data
	main:
		ch_data_expanded = ch_data.map {
			tuple ( it, it["normalBAM"], it["normalBAI"], it["tumorBAM"], it["tumorBAI"] )
		}
		delly_matched_call (genome_build, ch_fasta, ch_data_expanded)
		delly_matched_filter (genome_build, delly_matched_call.out.result)

	emit:
		results = delly_matched_filter.out.result
}

