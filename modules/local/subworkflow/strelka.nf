#!/usr/bin/env nextflow

include { interval_bed } from "../software/genome/main"
include { strelka_matched } from "../software/strelka/main"

workflow STRELKA {

	take:
		ch_fasta
		ch_dict
		ch_interval
		ch_data
		ch_indel

	main:
		ch_interval_list = ch_interval.collectFile (name: 'interval_names.tsv', newLine: true)
		interval_bed (ch_dict, ch_interval_list)

		ch_data_expanded = ch_data.map { it ->
			tuple ( it, it["normalBAM"], it["normalBAI"], it["tumorBAM"], it["tumorBAI"] )
		}

		strelka_matched (ch_fasta, interval_bed.out.result, ch_data_expanded, ch_indel)

	emit:
		result = strelka_matched.out.result

}

