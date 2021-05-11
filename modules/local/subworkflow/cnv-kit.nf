


include { cnv_kit_matched } from "../software/cnv-kit/main"


workflow CNV_KIT {
	take:
		ch_fasta
		ch_dict
		ch_interval
		ch_data

	main:
		ch_interval_list = ch_interval.collectFile (name: 'interval_names.tsv', newLine: true)
		interval_bed (ch_dict, ch_interval_list)

		ch_data_expanded = ch_data.map { it ->
			tuple ( it, it["normalBAM"], it["normalBAI"], it["tumorBAM"], it["tumorBAI"] )
		}

		cnv_kit_matched (ch_fasta, interval_bed.out.result, ch_data_expanded)

}

