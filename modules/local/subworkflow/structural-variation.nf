
include { structural_variation_matched } from "../software/structural-variation/main"

workflow STRUCTURAL_VARIATION
{

	take:
		ch_interval_bed
		ch_manta
		ch_delly

	main:
		ch_manta_keyed = ch_manta.map { tuple ( it[0]["sampleName"], it ) }
		ch_delly_keyed = ch_delly.map { tuple ( it[0]["sampleName"], it ) }

		ch_data = ch_manta_keyed.mix (ch_delly_keyed).groupTuple (size=2).view { "ch_data: ${it}" }


		structural_variation_matched (ch_interval_bed, ch_data)


}

