#!/usr/bin/env nextflow

include { bubble_tree_matched } from "../software/bubble-tree/main"


workflow BUBBLE_TREE {

	take:
		genome_build
		ch_interval_auto
		ch_ratio
		ch_loh

	main:
		ch_interval_csv_string = ch_interval_auto.map { it.join (",") }

		ch_loh_key = ch_loh.map { [it[0]["sampleName"], ["loh", [it[1]]], it[0]] }
		ch_ratio_key = ch_ratio.filter { it[1] == "Tumor" }
			.map { [it[0]["sampleName"], ["ratio", [it[2], it[4]]], it[0]] }

		ch_loh_and_ratio = ch_loh_key.mix (ch_ratio_key)
			.dump (tag: 'bubble tree before groupTuple')
			.groupTuple (size: 2)
			.map {
				def m = it[1].inject ([:]) { accumulator, item ->
					accumulator[item[0]] = item[1]
					accumulator
				}
				[it[2][0]] + m["loh"] + m["ratio"]
		}.dump (tag: 'bubble tree after groupTuple')

		bubble_tree_matched (genome_build, ch_interval_csv_string, ch_loh_and_ratio)

	emit:
		result = bubble_tree_matched.out.result
}


