#!/usr/bin/env nextflow

include { mutect_matched } from "../software/mutect/main"

workflow MUTECT
{
	take:
		ch_fasta
		ch_interval
		ch_data
	main:
		ch_data_interval = ch_data.map { it ->
			tuple ( it, it["normalBAM"], it["tumorBAM"] )
		}

		mutect_matched (ch_fasta, ch_interval, ch_data_interval)

	emit:
		results = mutect_matched.out.results
}

