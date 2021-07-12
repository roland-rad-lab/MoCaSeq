#!/usr/bin/env nextflow

include { bubble_tree_matched } from "../software/bubble-tree/main"


workflow BUBBLE_TREE {

	take:
		ch_ratio
		ch_loh

	main:
		ch_loh_key = ch_loh.map { [it[0]["sampleName"], ["loh", [it[1]]], it[0]] }
		ch_ratio_key = ch_ratio.filter { it[1] == "20000" }
			.map { [it[0]["sampleName"], ["ratio", [it[2], it[3]]], it[0]] }

		ch_loh_and_ratio = ch_loh_key.mix (ch_ratio_key)
			.groupTuple (size: 2)
			.map {
				def m = it[1].inject ([:]) { accumulator, item ->
					accumulator[item[0]] = item[1]
					accumulator
				}
				[it[2][0]] + m["loh"] + [m["ratio"][1]]
		}

		bubble_tree_matched (ch_loh_and_ratio)

	emit:
		result = bubble_tree_matched.out.result
}


