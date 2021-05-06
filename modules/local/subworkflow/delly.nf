#!/usr/bin/env nextflow

include { delly_matched; delly_matched_filter } from "../software/delly/main"

workflow DELLY
{
	take:
		ch_fasta
		ch_data
	main:
		delly_matched (ch_fasta, ch_data)
		delly_matched_filter (delly_matched.out.result)

	emit:
		results = delly_matched.out.result
}

