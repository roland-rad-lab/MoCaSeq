#!/usr/bin/env nextflow

params.hmm_copy = [:]

include { hmm_copy_wig as hmm_copy_wig_normal; hmm_copy_wig as hmm_copy_wig_tumor; hmm_copy_tsv; hmm_copy_plot } from "../software/hmm-copy/main"

workflow HMM_COPY {

	take:
		ch_interval
		ch_interval_bed
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
		ch_hmm_copy_wig_normal_keyed = hmm_copy_wig_normal.out.result.map { tuple ( tuple (it[0].sampleName, it[2]), it ) }
		ch_hmm_copy_wig_tumor_keyed = hmm_copy_wig_tumor.out.result.map { tuple ( tuple (it[0].sampleName, it[2]), it ) }

		ch_hmm_copy_wig = ch_hmm_copy_wig_normal_keyed.mix (ch_hmm_copy_wig_tumor_keyed).groupTuple (size:2)
			.map {
				def m = it[1].inject ([:]) { accumulator, item ->
					accumulator[item[1]] = [meta: item[0], resolution: item[2], wig: item[3]]
					accumulator
				}
				tuple ( m["Normal"]["resolution"], m["Normal"]["meta"], m["Normal"]["wig"], m["Tumor"]["wig"] )
			}
		ch_hmm_copy_wig_resolution = ch_wig_resolution.cross (ch_hmm_copy_wig)
			.dump (tag: 'hcw')
			.map {
				tuple ( it[0][0], it[0][1], it[0][2], it[1][1], it[1][2], it[1][3] )
			}

		hmm_copy_tsv (ch_interval_csv_string, ch_hmm_copy_wig_resolution)
		hmm_copy_plot (ch_interval_bed, hmm_copy_tsv.out.result)

	emit:
		tsv = hmm_copy_tsv.out.result
}



