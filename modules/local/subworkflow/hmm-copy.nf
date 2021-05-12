#!/usr/bin/env nextflow

params.hmm_copy = [:]

include { hmm_copy_wig as hmm_copy_wig_normal; hmm_copy_wig as hmm_copy_wig_tumor } from "../software/hmm-copy/main"

workflow HMM_COPY {

	take:
		ch_interval
		ch_gc_wig
		ch_map_wig
		ch_data

	main:
		ch_resolution = params.hmm_copy && params.hmm_copy.resolution && params.hmm_copy.resolution ? params.hmm_copy.resolution.tokenize (",") : Channel.empty ()

		ch_interval_csv_string = ch_interval.toList ().map { it.join (",") }
		ch_data_expanded_normal = ch_data.map { tuple (it, "Normal", it["normalBAM"], it["normalBAI"] ) }
		ch_data_expanded_tumor = ch_data.map { tuple (it, "Tumor", it["tumorBAM"], it["tumorBAI"] ) }

		hmm_copy_wig_normal (ch_interval_csv_string, ch_resolution, ch_data_expanded_normal)
		hmm_copy_wig_tumor (ch_interval_csv_string, ch_resolution, ch_data_expanded_tumor)

		ch_gc_wig_resolution = ch_gc_wig.map { tuple ( it.name , it) }.view ()
}


