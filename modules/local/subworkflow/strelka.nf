#!/usr/bin/env nextflow

include { strelka_matched; strelka_matched_post } from "../software/strelka/main"

workflow STRELKA {

	take:
		genome_build
		ch_fasta
		ch_interval_bed
		ch_data
		ch_indel

	main:
		ch_indel_branched = ch_indel.branch {
			matched: it[1] == "matched"
			other: true
		}

		ch_indel_branched.other.view { "[MoCaSeq] error: Failed to find matching STRELKA workflow path for input:\n${it}" }

		ch_data_expanded_matched = ch_indel_branched.matched.map { it ->
			tuple ( it[0], it[0]["normalBAM"], it[0]["normalBAI"], it[0]["tumorBAM"], it[0]["tumorBAI"], it[2], it[3] )
		}

		strelka_matched (genome_build, ch_fasta, ch_interval_bed, ch_data_expanded_matched)
		strelka_matched_post (genome_build,strelka_matched.out.result)

	emit:
		result = strelka_matched.out.result

}

