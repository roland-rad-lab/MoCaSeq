#!/usr/bin/env nextflow

include { igv_track_depth as igv_track_depth_normal } from "../software/igv-track/main"

workflow IGV_TRACK_READ {

	take:
		ch_interval
		ch_interval_bed
		ch_data

	main:
		ch_interval_space_string = ch_interval.toList ().map { it.join (" ") }
		ch_data_normal = ch_data.map { tuple (it, it["normalBAM"], it["normalBAI"] ) }
		ch_data_tumor = ch_data.map { tuple (it, it["tumorBAM"], it["tumorBAI"] ) }

		igv_track_depth_normal (ch_interval_space_string, ch_interval_bed, ch_data_normal)
}




