#!/usr/bin/env nextflow

include {
	igv_track_depth as igv_track_depth_normal;
	igv_track_depth as igv_track_depth_tumor;
	igv_track_cn
} from "../software/igv-track/main"

workflow IGV_TRACK_READ {

	take:
		ch_interval
		ch_interval_bed
		ch_data

	main:
		ch_interval_space_string = ch_interval.toList ().map { it.join (" ") }
		ch_data_normal = ch_data.map { tuple (it, "Normal", it["normalBAM"], it["normalBAI"] ) }
		ch_data_tumor = ch_data.map { tuple (it, "Tumor", it["tumorBAM"], it["tumorBAI"] ) }

		igv_track_depth_normal (ch_interval_space_string, ch_interval_bed, ch_data_normal)
		igv_track_depth_tumor (ch_interval_space_string, ch_interval_bed, ch_data_tumor)
}

workflow IGV_TRACK_CN {

	take:
		ch_cns

	main:
		igv_track_cn (ch_cns)
}


