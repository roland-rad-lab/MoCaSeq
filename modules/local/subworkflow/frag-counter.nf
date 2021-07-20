#!/usr/bin/env nextflow

params.frag_counter = [:]

workflow FRAG_COUNTER {

	take:
		ch_interval
		ch_gc_wig
		ch_map_wig
		ch_data

	main:
		ch_resolution = params.hmm_copy && params.hmm_copy.resolution && params.hmm_copy.resolution ? Channel.value (params.frag_counter.resolution.tokenize (",")) : Channel.empty ()

		ch_interval_csv_string = ch_interval.toList ().map { it.join (",") }
		ch_data_expanded_normal = ch_data.map { tuple (it, "Normal", it["normalBAM"], it["normalBAI"] ) }

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
			.dump (tag: 'frag_counter ch_wig_resolution')

}

