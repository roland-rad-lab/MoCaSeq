#!/usr/bin/env nextflow

include { interval_bed } from "../software/genome/main"
include {
	bash_expand_path as bash_expand_path_gc
	bash_expand_path as bash_expand_path_map
} from "../helpers"

params.genomes = [:]
params.genome_annotations = [:]

workflow PREPARE_GENOME
{

	take:
		genome_name

	main:
		if ( genome_name == null ) { exit 1, "[MoCaSeq] error: Genome name not found. Check params.genome_build." }
		ch_fasta = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["fasta"] ? Channel.of (params.genomes[genome_name]["fasta"]).first () : Channel.empty ()
		ch_fasta_index = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["fasta_index"] ? Channel.of (params.genomes[genome_name]["fasta_index"]).first () : Channel.empty ()
		ch_fasta_index_flat = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["fasta_index_flat"] ? Channel.of (params.genomes[genome_name]["fasta_index_flat"]).first () : Channel.empty ()
		ch_dict = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["dict"] ? Channel.of (params.genomes[genome_name]["dict"]).first () : Channel.empty ()
		ch_dir = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["dir"] ? Channel.of (params.genomes[genome_name]["dir"]).first () : Channel.empty ()
		ch_chrom_names = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["names"] && params.genomes[genome_name]["names"]["auto_sex"] ? Channel.fromList (params.genomes[genome_name]["names"]["auto_sex"]) : Channel.empty ()
		chrom_n = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["names"] && params.genomes[genome_name]["names"]["auto_sex"] ? params.genomes[genome_name]["names"]["auto_sex"].size () : 0

		ch_interval_list = ch_chrom_names.collectFile (name: 'interval_names.tsv', newLine: true)
		interval_bed (ch_dict, ch_interval_list)

	emit:
		fasta            = ch_fasta
		fasta_index      = ch_fasta_index
		fasta_index_flat = ch_fasta_index_flat
		dict             = ch_dict
		dir              = ch_dir
		chrom_names      = ch_chrom_names
		interval_bed     = interval_bed.out.result.first ()
		_chrom_n         = chrom_n
		
}

workflow GENOME_ANNOTATION
{
	take:
		genome_name

	main:
		if ( genome_name == null ) { exit 1, "[MoCaSeq] error: Genome name not found. Check params.genome_build." }
		ch_gc_wig = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["gc_wig"] ? Channel.of (params.genome_annotations[genome_name]["gc_wig"]) : Channel.empty ()
		ch_map_wig = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["map_wig"] ? Channel.of (params.genome_annotations[genome_name]["map_wig"]) : Channel.empty ()
		ch_common_vcf = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["common_vcf"] ? Channel.of (params.genome_annotations[genome_name]["common_vcf"]).first () : Channel.empty ()
		ch_gencode_genes_bed = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["gencode_genes_bed"] ? Channel.of (params.genome_annotations[genome_name]["gencode_genes_bed"]).first () : Channel.empty ()
		ch_micro_satellite = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["micro_satellite"] ? Channel.of (params.genome_annotations[genome_name]["micro_satellite"]).first () : Channel.empty ()

		bash_expand_path_gc (ch_gc_wig)
		bash_expand_path_map (ch_map_wig)

	emit:
		gc_wig = bash_expand_path_gc.out.splitText ()
		map_wig = bash_expand_path_map.out.splitText ()
		common_vcf = ch_common_vcf
		gencode_genes_bed = ch_gencode_genes_bed
		micro_satellite = ch_micro_satellite
}


