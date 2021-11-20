#!/usr/bin/env nextflow

params.hmm_copy = [:]

include { hmm_copy_wig as hmm_copy_wig_normal; hmm_copy_wig as hmm_copy_wig_tumor; hmm_copy_tsv; hmm_copy_plot } from "../software/hmm-copy/main"

workflow HMM_COPY {

	take:
		genome_build
		ch_interval
		ch_interval_bed
		ch_gc_wig
		ch_map_wig
		ch_data

	main:
		ch_resolution = params.hmm_copy && params.hmm_copy.resolution && params.hmm_copy.resolution ? params.hmm_copy.resolution.tokenize (",") : Channel.empty ()
		ch_interval_csv_string = ch_interval.map { it.join (",") }

		ch_data_branched = ch_data.branch {
			normal: it["type"] == "Normal"
			tumor: it["type"] == "Tumor"
			other: true
		}

		ch_data_branched.other.view { "[MoCaSeq] error: Unknown (type) for input:\n${it}\nExpected: [Normal,Tumor]." }

		ch_data_expanded_normal = ch_data_branched.normal.map { tuple (it, "Normal", it["normalBAM"], it["normalBAI"] ) }
		ch_data_expanded_tumor = ch_data_branched.tumor.map { tuple (it, "Tumor", it["tumorBAM"], it["tumorBAI"] ) }

		hmm_copy_wig_normal (genome_build, ch_interval_csv_string, ch_resolution, ch_data_expanded_normal)
		hmm_copy_wig_tumor (genome_build, ch_interval_csv_string, ch_resolution, ch_data_expanded_tumor)

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
		ch_hmm_copy_wig_normal_keyed = hmm_copy_wig_normal.out.result.map { tuple ( tuple (groupKey (it[0].sampleGroup, it[0].sampleGroupSize), it[2]), it ) }
		ch_hmm_copy_wig_tumor_keyed = hmm_copy_wig_tumor.out.result.map { tuple ( tuple (groupKey (it[0].sampleGroup, it[0].sampleGroupSize), it[2]), it ) }

		ch_hmm_copy_wig = ch_hmm_copy_wig_normal_keyed.mix (ch_hmm_copy_wig_tumor_keyed).groupTuple (remainder: true)
			.flatMap {
				// Collect [Normal,Tumor] samples in the Sample_Group, output tuple (with normal) for each Tumor sample
				def m = it[1].inject ([:]) { accumulator, item ->
					if ( accumulator.containsKey (item[1]) )
					{
						accumulator[item[1]].add ([meta: item[0], resolution: item[2], wig: item[3]])
					}
					else
					{
						accumulator[item[1]] = [[meta: item[0], resolution: item[2], wig: item[3]]]
					}
					accumulator
				}
				if ( m.containsKey ("Normal") && m.containsKey ("Tumor") )
				{
					m["Tumor"].collect { jt -> jt.put ("wigNormal", m["Normal"][0]["wig"]);[jt["resolution"], jt] }
				}
				else
				{
					[]
				}
			}
		ch_hmm_copy_wig_resolution = ch_wig_resolution.cross (ch_hmm_copy_wig)
			.map {
				tuple ( it[0][0], it[0][1], it[0][2], it[1][1]["meta"], it[1][1]["wig"], it[1][1]["wigNormal"] )
			}.dump (tag: 'hmm-copy wr')

		hmm_copy_tsv (genome_build, ch_interval_csv_string, ch_hmm_copy_wig_resolution)
		hmm_copy_plot (genome_build, ch_interval_bed, hmm_copy_tsv.out.result)

	emit:
		cnr = hmm_copy_tsv.out.cnr
		call = hmm_copy_tsv.out.call
}



