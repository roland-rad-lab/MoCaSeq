#!/usr/bin/env nextflow

include { interval_bed_intersect } from "../software/genome/main"
include { cnv_kit_matched; cnv_kit_single; cnv_kit_segment; cnv_kit_target_bed; cnv_kit_coverage; cnv_kit_fix; cnv_kit_reference } from "../software/cnv-kit/main"

workflow CNV_KIT {

	take:
		genome_build
		ch_fasta
		ch_fasta_index_flat
		ch_interval_bed
		ch_gencode_genes_bed
		ch_data

	main:
		interval_bed_intersect (ch_gencode_genes_bed, ch_interval_bed.map { it[0] }, Channel.of ("-wa"))

		ch_interval_bed_intersection = interval_bed_intersect.out.result.first ()

		ch_data_expanded = ch_data.filter { it["type"] == "Tumor" }.map { it ->
			tuple (it, it["normalBAM"], it["normalBAI"], it["tumorBAM"], it["tumorBAI"] )
		}

		ch_data_branched = ch_data.branch {
			normal: it["type"] == "Normal"
			tumor: it["type"] == "Tumor"
			other: true
		}

		ch_data_branched.other.view { "[MoCaSeq] error: Unknown (type) for input:\n${it}\nExpected: [Normal,Tumor]." }

		ch_data_expanded_normal = ch_data_branched.normal.map { it -> tuple (it, "Normal", it["normalBAM"], it["normalBAI"]) }
		ch_data_expanded_tumor = ch_data_branched.tumor.map { it -> tuple (it, "Tumor", it["tumorBAM"], it["tumorBAI"]) }

		cnv_kit_matched (genome_build, ch_fasta, ch_fasta_index_flat, ch_interval_bed_intersection, ch_data_expanded)
		cnv_kit_single (genome_build, ch_fasta, ch_fasta_index_flat, ch_interval_bed_intersection, ch_data_expanded_normal.mix (ch_data_expanded_tumor))

	emit:
		cns = cnv_kit_single.out.cns
}

workflow CNV_KIT_FIX {

	take:
		genome_build
		ch_reference
		ch_coverage

	main:
		cnv_kit_fix (genome_build, ch_reference, ch_coverage)

	emit:
		result = cnv_kit_fix.out.result
}

workflow CNV_KIT_SEGMENT {

	take:
		genome_build
		ch_coverage
	main:
		cnv_kit_segment (genome_build, ch_coverage)

	emit:
		cns = cnv_kit_segment.out.cns
		call = cnv_kit_segment.out.call
}

workflow CNV_KIT_SELECT_SAMPLE {


	take:
		genome_build
		size_tsv

	main:
		ch_representative_sample = size_tsv.splitCsv (header: true, sep: "\t")
		.filter {
			if ( !it.containsKey ("sample") ) exit 1, "[MoCaSeq] error: Invalid lines in normal_size_tsv, no 'sample' column '${it}'"
			if ( !it.containsKey ("genome_build") ) exit 1, "[MoCaSeq] error: Invalid lines in normal_size_tsv, no 'genome_build' column '${it}'"
			if ( !it.containsKey ("bam_size") ) exit 1, "[MoCaSeq] error: Invalid lines in normal_size_tsv, no 'bam_size' column '${it}'"
			if ( !it["bam_size"].isLong () ) exit 1, "[MoCaSeq] error: Invalid lines in normal_size_tsv, 'bam_size' column was not a long integer '${it}'"
			return it["genome_build"] == genome_build
		}.toSortedList( { a, b -> b["bam_size"] as long <=> a["bam_size"] as long } )
		.map {
			if ( it.isEmpty () ) exit 1, "[MoCaSeq] error: No samples were found in normal_size_tsv"
			["sample", "genome_build"].collect { jt -> it[it.size ().intdiv (2)][jt] }
		}

	emit:
		result = ch_representative_sample
}

workflow CNV_KIT_TARGET_BED {

	take:
		genome_build
		ch_fasta
		ch_interval
		ch_data_expanded

	main:
		ch_interval_csv_string = ch_interval.map { it.join (",") }
		cnv_kit_target_bed (genome_build, ch_fasta, ch_interval_csv_string, ch_data_expanded)

		ch_parsed = cnv_kit_target_bed.out.result.map {
			def jsonSlurper = new groovy.json.JsonSlurper ()
			def data = jsonSlurper.parse (it[1])
			tuple ( it[0], data["target_avg_size"], data["wgs_depth"] )
		}
	emit:
		result = ch_parsed.first ()
}

workflow CNV_KIT_COVERAGE {

	take:
		genome_build
		ch_fasta
		ch_interval_bed
		ch_data_expanded
	main:
		cnv_kit_coverage (genome_build, ch_fasta, ch_interval_bed, ch_data_expanded)

	emit:
		result = cnv_kit_coverage.out.result
}

workflow CNV_KIT_PON {

	take:
		genome_build
		ch_fasta
		ch_normal_coverage_tsv

	main:

		ch_normal_coverage_tsv_filtered = ch_normal_coverage_tsv.splitCsv (header: true, sep: "\t")
		.filter {
			if ( !it.containsKey ("genome_build") ) exit 1, "[MoCaSeq] error: Invalid lines in normal_coverage_tsv, no 'genome_build' column '${it}'"
			return it["genome_build"] == genome_build
		}.tap { ch_paths }
		.map { jt ->
			def f = file (jt["normal_cov"],glob: false)
			// DO NOT MODIFY jt
			def header = ["sample", "genome_build", "normal_cov"]
			def ldata = [jt["sample"], jt["genome_build"], f.name]
			[ header.join ("\t"), ldata.join ("\t"), ""].join ("\n")
		}
		.collectFile (keepHeader: true, skip: 1)

		cnv_kit_reference (genome_build, ch_fasta, ch_paths.map { file (it["normal_cov"], glob: false) }.collect (), ch_normal_coverage_tsv_filtered)
}

