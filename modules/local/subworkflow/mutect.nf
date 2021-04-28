#!/usr/bin/env nextflow

include { mutect_matched } from "../software/mutect/main"

workflow MUTECT
{
	take:
		genome
		data
	main:
		ch_data_chrom = data.map { it ->
			tuple ( it, it["NormalBAM"], it["TumorBAM"] )
		}.combine (genome["chrom_names"])

		mutect_matched (ch_data_chrom)

	emit:
		results = mutect_matched.out.results
}

