#!/usr/bin/env nextflow

include { interval_bed; interval_bed_intersect } from "../software/genome/main"
include { cnv_kit_matched; cnv_kit_single as cnv_kit_single_normal; cnv_kit_single as cnv_kit_single_tumor } from "../software/cnv-kit/main"


workflow CNV_KIT {
	take:
		ch_fasta
		ch_fasta_index_flat
		ch_dict
		ch_interval
		ch_gencode_genes_bed
		ch_data

	main:
		ch_interval_list = ch_interval.collectFile (name: 'interval_names.tsv', newLine: true)
		interval_bed (ch_dict, ch_interval_list)
		interval_bed_intersect (ch_gencode_genes_bed, interval_bed.out.result.map { it[0] }, Channel.of ("-wa"))

		ch_data_expanded = ch_data.map { it ->
			tuple (it, it["normalBAM"], it["normalBAI"], it["tumorBAM"], it["tumorBAI"] )
		}
		ch_data_expanded_normal = ch_data.map { it ->
			tuple (it, "Normal", it["normalBAM"], it["normalBAI"])
		}
		ch_data_expanded_tumor = ch_data.map { it ->
			tuple (it, "Tumor", it["tumorBAM"], it["tumorBAI"])
		}

		cnv_kit_matched (ch_fasta, ch_fasta_index_flat, interval_bed_intersect.out.result, ch_data_expanded)
		cnv_kit_single_normal (ch_fasta, ch_fasta_index_flat, interval_bed_intersect.out.result, ch_data_expanded_normal)
		cnv_kit_single_tumor (ch_fasta, ch_fasta_index_flat, interval_bed_intersect.out.result, ch_data_expanded_tumor)
}

