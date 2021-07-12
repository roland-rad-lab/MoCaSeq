
include { jabba_matched;jabba_plot } from "../software/jabba/main"

workflow JABBA
{
	take:
		ch_manta
		ch_ratio

	main:

		ch_manta_key = ch_manta.map { [it[0]["sampleName"], ["manta", [it[1]]], it[0]] }
		ch_ratio_key = ch_ratio.filter { it[1] == "1000" }
			.map { [it[0]["sampleName"], ["ratio", [it[2], it[3]]], it[0]] }

		ch_manta_and_ratio = ch_manta_key.mix (ch_ratio_key)
			.groupTuple (size: 2)
			.map {
				def m = it[1].inject ([:]) { accumulator, item ->
					accumulator[item[0]] = item[1]
					accumulator
				}
				[it[2][0]] + m["manta"] + m["ratio"]
		}

		jabba_matched (ch_manta_and_ratio)
		jabba_plot (jabba_matched.out.result)
}


