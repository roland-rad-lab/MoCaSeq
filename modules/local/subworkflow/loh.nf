#!/usr/bin/env nextflow

include { loh_matched; loh_matched_assign_alleles; } from "../software/loh/main"

workflow LOH {

	take:
		genome_build
		ch_fasta
		ch_fasta_index
		ch_interval
		ch_interval_bed
		ch_data

	main:
		ch_interval_csv_string = ch_interval.map { it.join (",") }

		ch_data_branched = ch_data.branch {
			single: it[1] == "Normal" || it[1] == "Tumor"
			matched: it[1] == "matched"
			other: true
		}

		ch_data_branched.other.view { "[MoCaSeq] error: Failed to find matching LOH workflow path for input:\n${it}" }

		ch_data_single_sample = ch_data_branched.single.map { [groupKey (it[0].sampleGroup, it[0].sampleGroupSize), it] }
			.groupTuple (remainder: true)
			.flatMap {
				// Collect [Normal,Tumor] samples in the Sample_Group, output tuple (with normal) for each Tumor sample
				def m = it[1].inject ([:]) { accumulator, item ->
					if ( accumulator.containsKey (item[1]) )
					{
						accumulator[item[1]].add ([meta: item[0], vcf: item[2], vcf_index: item[3]])
					}
					else
					{
						accumulator[item[1]] = [[meta: item[0], vcf: item[2], vcf_index: item[3]]]
					}
					accumulator
				}
				if ( m.containsKey ("Normal") && m.containsKey ("Tumor") )
				{
					m["Tumor"].collect { jt -> [jt["meta"], m["Normal"][0]["vcf"], m["Normal"][0]["vcf_index"], jt["vcf"], jt["vcf_index"]] }
				}
				else
				{
					[]
				}
			}

		loh_matched (genome_build, ch_interval_bed, ch_data_single_sample)
		loh_matched_assign_alleles (genome_build, ch_fasta, ch_fasta_index, ch_interval_csv_string, loh_matched.out.result)

	emit:
		result = loh_matched_assign_alleles.out.result
}

