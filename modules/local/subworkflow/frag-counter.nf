#!/usr/bin/env nextflow

params.frag_counter = [:]

include { frag_counter_wig_to_rds;frag_counter } from "../software/frag-counter/main"

workflow FRAG_COUNTER {

	take:
		genome_name
		ch_interval
		ch_gc_wig
		ch_map_wig
		ch_data

	main:
		ch_resolution = params.frag_counter && params.frag_counter["resolution"] ? Channel.value (params.frag_counter["resolution"].toString ()) : Channel.empty ()

		ch_interval_csv_string = ch_interval.toList ().map { it.join (",") }

		ch_gc_wig_resolution = ch_gc_wig.map { new File (it.trim ()) }.map {
			def m = (it.name =~ /\.([\w\-]+)\.wig$/)
			if ( m.count != 1 ) exit 1, "[MoCaSeq] error: Failed to parse resolution from wig file '${it.name}'"
			tuple (m[0][1] , it.path)
		}
		ch_map_wig_resolution = ch_map_wig.map { new File (it.trim ()) }.map {
			def m = (it.name =~ /\.([\w\-]+)\.wig$/)
			if ( m.count != 1 ) exit 1, "[MoCaSeq] error: Failed to parse resolution from wig file '${it.name}'"
			tuple (m[0][1], it.path)
		}

		ch_wig_resolution = ch_gc_wig_resolution.join (ch_map_wig_resolution)
			.join (ch_resolution.flatten ().map { [it] })

		ch_data_expanded_normal = ch_data.map { it ->
			tuple (it, "Normal", it["normalBAM"], it["normalBAI"])
		}
		ch_data_expanded_tumor = ch_data.map { it ->
			tuple (it, "Tumor", it["tumorBAM"], it["tumorBAI"])
		}

		ch_data_expanded = ch_data_expanded_normal.mix (ch_data_expanded_tumor)

		frag_counter_wig_to_rds (genome_name, ch_interval_csv_string, ch_wig_resolution)
		frag_counter (frag_counter_wig_to_rds.out.result.first (), ch_data_expanded)

	emit:
		result = frag_counter.out.result
}
