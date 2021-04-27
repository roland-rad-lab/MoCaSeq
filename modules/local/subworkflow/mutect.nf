#!/usr/bin/env nextflow

include {  } from "../software/mutect/main.nf"

workflow mutect
{
	take:
		genome
		data
	main:
		ch_data_chrom = data.map { it ->
			tuple ( it, it["NormalBAM"], it["TumorBAM"] )
		}.combine (genome.chrom_names.auto_sex)

		mutect_matched (ch_data_chrom)

	emit:
		mutect_matched.out
}

