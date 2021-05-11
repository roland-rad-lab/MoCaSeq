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

		ch_indel_branched = ch_indel.branch {
			matched: it[1] == "matched"
			other: true
		}

		ch_indel_branched.other.view { "[MoCaSeq] error: Failed to find matching STRELKA workflow path for input:\n${it}" }

		ch_data_expanded_matched = ch_indel_branched.matched.map { it ->
			tuple ( it, it["normalBAM"], it["normalBAI"], it["tumorBAM"], it["tumorBAI"], it[2], it[3] )
		}

		strelka_matched (ch_fasta, interval_bed.out.result, ch_data_expanded_matched)

	emit:
		result = strelka_matched.out.result

}

