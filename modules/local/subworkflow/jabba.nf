
include { jabba_matched;jabba_plot } from "../software/jabba/main"

workflow JABBA
{
	take:
		genome_build
		ch_interval
		ch_manta
		ch_coverage
		ch_segment
		ch_bubble

	main:

		ch_interval_csv_string = ch_interval.map { it.join (",") }
		ch_manta_key = ch_manta.map { [it[0]["sampleName"], ["manta", [it[1]]], it[0]] }
		ch_coverage_key = ch_coverage.filter { it[1] == "Tumor" }
			.map { [it[0]["sampleName"], ["coverage", [it[2], it[4]]], it[0]] }
		ch_segment_key = ch_segment.filter { it[1] == "Tumor" }
			.map { [it[0]["sampleName"], ["segment", [it[2], it[4]]], it[0]] }
		ch_bubble_key = ch_bubble.map {
				def bubble_tree_output = it[1].getText ()
				def m_bubble_tree_output = bubble_tree_output =~ /\s*(\w+):\s+?([0-9\.,\s]+);?/
				def result = m_bubble_tree_output.iterator ().toList ().inject ([:]) { accumulator, item ->
					accumulator[item[1]] = item[2]
					accumulator
				}
				[it[0]["sampleName"], ["bubble", [result["Purity"], result["Ploidy"]]], it[0]]
			}

		ch_manta_and_ratio_and_bubble = ch_manta_key.mix (ch_coverage_key,ch_segment_key,ch_bubble_key)
			.dump (tag: 'jabba before groupTuple')
			.groupTuple (size: 4)
			.map {
				def m = it[1].inject ([:]) { accumulator, item ->
					accumulator[item[0]] = item[1]
					accumulator
				}
				[it[2][0]] + m["manta"] + m["coverage"] + m["segment"] + m["bubble"]
			}
			.dump (tag: 'manta_ratio_and_bubble')

		jabba_matched (genome_build, ch_manta_and_ratio_and_bubble)
		jabba_plot (genome_build, ch_interval_csv_string, jabba_matched.out.result)
}


