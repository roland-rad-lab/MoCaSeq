
include { structural_variation_matched } from "../software/structural-variation/main"

workflow STRUCTURAL_VARIATION
{

	take:
		ch_interval_bed
		ch_manta
		ch_delly
		ch_cnv_kit

	main:
		ch_manta_keyed = ch_manta.map { tuple ( it[0]["sampleName"], "Manta", it ) }
		ch_delly_keyed = ch_delly.map { tuple ( it[0]["sampleName"], "Delly", it ) }
		ch_cnv_kit_keyed = ch_cnv_kit.map { tuple ( it[0]["sampleName"], "CNVKit", it ) }

		ch_data = ch_manta_keyed.mix (ch_delly_keyed,ch_cnv_kit_keyed).groupTuple ().view { "ch_data: ${it}" }
			.map {
				def manta_index = it[1].indexOf ("Manta")
				def delly_index = it[1].indexOf ("Delly")
				def cnvkit_index = it[1].indexOf ("CNVKit")

				def manta_vcf = it[2][manta_index][1]
				def manta_vcf_index = it[2][manta_index][2]
				def delly_bcf = it[2][delly_index][1]
				def cnvkit_normal = it[2][cnvkit_index][1]
				def cnvkit_tumor = it[2][cnvkit_index][2]

				tuple ( ["sampleName": it[0]], manta_vcf, manta_vcf_index, delly_bcf, cnvkit_normal, cnvkit_tumor)
			}.view { "ch_data out: ${it}" }


		structural_variation_matched (ch_interval_bed, ch_data)

}

