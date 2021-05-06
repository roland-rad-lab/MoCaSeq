#!/usr/bin/env nextflow

include { delly_matched_call; delly_matched_filter } from "../software/delly/main"

workflow DELLY
{
	take:
		ch_fasta
		ch_data
	main:
		delly_matched_call (ch_fasta, ch_data)
		delly_matched_filter (delly_matched_call.out.result)

	emit:
		results = delly_matched_filter.out.result
}

