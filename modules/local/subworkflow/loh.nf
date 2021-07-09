#!/usr/bin/env nextflow

include { loh_matched; loh_matched_assign_alleles; } from "../software/loh/main"

workflow LOH {

	take:
		ch_fasta
		ch_fasta_index
		ch_interval_bed
		ch_data

	main:
		ch_data_branched = ch_data.branch {
			single: it[1] == "Normal" || it[1] == "Tumor"
			matched: it[1] == "matched"
			other: true
		}

		ch_data_branched.other.view { "[MoCaSeq] error: Failed to find matching LOH workflow path for input:\n${it}" }

		ch_data_single_sample = ch_data_branched.single.map { [it[0]["sampleName"], it] }
			.groupTuple (size: 2)
			.map { it[1] }
			.dump (tag: 'LOH after groupTuple')
			.map {
				def m = it.inject ([:]) { accumulator, item ->
					accumulator[item[1]] = [item[2],item[3]]
					accumulator
				}
				[it[0][0]] + m["Normal"] + m["Tumor"]
			}
			.dump (tag: 'LOH after map')

		loh_matched (ch_fasta, ch_fasta_index, ch_interval_bed, ch_data_single_sample)
		loh_matched_assign_alleles (loh_matched.out.result)
}

