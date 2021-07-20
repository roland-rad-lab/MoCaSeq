#!/usr/bin/env nextflow

include { msi_matched } from "../software/msi-sensor/main"


workflow MSI_SENSOR {

	take:
		ch_micro_satellite
		ch_data

	main:
		ch_data_expanded = ch_data.map { it ->
			tuple ( it, it["normalBAM"], it["normalBAI"], it["tumorBAM"], it["tumorBAI"] )
		}

		msi_matched (ch_micro_satellite, ch_data_expanded)

	emit:
		result = msi_matched.out.result
}

