#!/usr/bin/env nextflow

include { interval_bed; manta_matched } from "../software/manta/main"

workflow MANTA
{
	take:
		ch_fasta
		ch_dict
		ch_interval
		ch_data

	main:
		ch_interval_list = ch_interval.collectFile (name: 'interval_names.tsv', newLine: true)
		interval_bed (ch_dict, ch_interval_list)

		ch_data_expanded = ch_data.map { it ->
			tuple ( it, it["normalBAM"], it["tumorBAM"] )
		}

		manta_matched (ch_fasta, interval_bed.out.result, ch_data_expanded)

	emit:
		result = manta_matched.out.result		
}

