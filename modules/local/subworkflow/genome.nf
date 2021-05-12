#!/usr/bin/env nextflow

params.genomes = [:]

workflow PREPARE_GENOME
{

	take:
		genome_name

	main:
		ch_fasta = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["fasta"] ? Channel.fromPath (params.genomes[genome_name]["fasta"]) : Channel.empty ()
		ch_fasta_index_flat = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["fasta_index_flat"] ? Channel.fromPath (params.genomes[genome_name]["fasta_index_flat"]) : Channel.empty ()
		ch_dict = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["dict"] ? Channel.fromPath (params.genomes[genome_name]["dict"]) : Channel.empty ()
		ch_chrom_names = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["names"] && params.genomes[genome_name]["names"]["auto_sex"] ? Channel.fromList (params.genomes[genome_name]["names"]["auto_sex"]) : Channel.empty ()
		chrom_n = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["names"] && params.genomes[genome_name]["names"]["auto_sex"] ? params.genomes[genome_name]["names"]["auto_sex"].size () : 0

	emit:
		fasta            = ch_fasta
		fasta_index_flat = ch_fasta_index_flat
		dict             = ch_dict
		chrom_names      = ch_chrom_names
		_chrom_n         = chrom_n
		
}

workflow GENOME_ANNOTATION
{
	take:
		genome_name

	main:
		ch_gencode_genes_bed = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["gencode_genes_bed"] ? Channel.fromPath (params.genome_annotations[genome_name]["gencode_genes_bed"]) : Channel.empty ()

	emit:
		gencode_genes_bed = ch_gencode_genes_bed
}


