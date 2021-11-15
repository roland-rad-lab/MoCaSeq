#!/usr/bin/env nextflow

include { msi_matched } from "../software/msi-sensor/main"


workflow MSI_SENSOR {

	take:
		genome_build
		ch_micro_satellite
		ch_data

	main:
		ch_data_expanded = ch_data.filter { it["type"] == "Tumor" }.map { it ->
			tuple ( it, it["normalBAM"], it["normalBAI"], it["tumorBAM"], it["tumorBAI"] )
		}

		msi_matched (genome_build, ch_micro_satellite, ch_data_expanded)

	emit:
		result = msi_matched.out.result
}

