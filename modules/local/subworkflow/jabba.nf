
include { jabba_matched;jabba_plot } from "../software/jabba/main"

workflow JABBA
{
	take:
		ch_manta
		ch_ratio
		ch_bubble

	main:

		ch_manta_key = ch_manta.map { [it[0]["sampleName"], ["manta", [it[1]]], it[0]] }
		ch_ratio_key = ch_ratio.filter { it[1] == "Tumor" && it[3] == "1000" }
			.map { [it[0]["sampleName"], ["ratio", [it[2], it[4], it[5]]], it[0]] }
		ch_bubble_key = ch_bubble.map {
				def bubble_tree_output = it[1].getText ()
				def m_bubble_tree_output = bubble_tree_output =~ /\s*(\w+):\s+?([0-9\.,\s]+);?/
				def result = m_bubble_tree_output.iterator ().toList ().inject ([:]) { accumulator, item ->
					accumulator[item[1]] = item[2]
					accumulator
				}
				[it[0]["sampleName"], ["bubble", [result["Purity"], result["Ploidy"]]], it[0]]
			}

		ch_manta_and_ratio_and_bubble = ch_manta_key.mix (ch_ratio_key,ch_bubble_key)
			.dump (tag: 'jabba before groupTuple')
			.groupTuple (size: 3)
			.map {
				def m = it[1].inject ([:]) { accumulator, item ->
					accumulator[item[0]] = item[1]
					accumulator
				}
				[it[2][0]] + m["manta"] + m["ratio"] + m["bubble"]
			}
			.dump (tag: 'manta_ratio_and_bubble')

		jabba_matched (ch_manta_and_ratio_and_bubble)
		jabba_plot (jabba_matched.out.result)
}


