#!/usr/bin/env nextflow

include { interval_bed_intersect } from "../software/genome/main"
include { cnv_kit_matched; cnv_kit_single as cnv_kit_single_normal; cnv_kit_single as cnv_kit_single_tumor; cnv_kit_segment } from "../software/cnv-kit/main"

workflow CNV_KIT {

	take:
		ch_fasta
		ch_fasta_index_flat
		ch_interval_bed
		ch_gencode_genes_bed
		ch_data

	main:
		interval_bed_intersect (ch_gencode_genes_bed, ch_interval_bed.map { it[0] }, Channel.of ("-wa"))

		ch_interval_bed_intersection = interval_bed_intersect.out.result.first ()

		ch_data_expanded = ch_data.map { it ->
			tuple (it, it["normalBAM"], it["normalBAI"], it["tumorBAM"], it["tumorBAI"] )
		}
		ch_data_expanded_normal = ch_data.map { it ->
			tuple (it, "Normal", it["normalBAM"], it["normalBAI"])
		}
		ch_data_expanded_tumor = ch_data.map { it ->
			tuple (it, "Tumor", it["tumorBAM"], it["tumorBAI"])
		}

		cnv_kit_matched (ch_fasta, ch_fasta_index_flat, ch_interval_bed_intersection, ch_data_expanded)
		cnv_kit_single_normal (ch_fasta, ch_fasta_index_flat, ch_interval_bed_intersection, ch_data_expanded_normal)
		cnv_kit_single_tumor (ch_fasta, ch_fasta_index_flat, ch_interval_bed_intersection, ch_data_expanded_tumor)

	emit:
		cns_normal = cnv_kit_single_normal.out.cns
		cns_tumor = cnv_kit_single_tumor.out.cns
}

workflow CNV_KIT_SEGMENT {

	take:
		coverage_source
		ch_coverage

	main:
		cnv_kit_segment (coverage_source, ch_coverage.tap { ch_coverage_copy })

	emit:
		tsv = cnv_kit_segment.out.result
}

