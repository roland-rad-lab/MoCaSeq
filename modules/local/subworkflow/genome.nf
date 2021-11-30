#!/usr/bin/env nextflow

include {
	interval_bed
	cache_genome_url as cache_genome_url_bwa_index
	cache_genome_url as cache_genome_url_dict
	cache_genome_url as cache_genome_url_fasta
	cache_genome_url as cache_genome_url_fasta_index
	cache_genome_url as cache_genome_url_common_vcf
} from "../software/genome/main"
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
		ch_ext_bwa_index = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["ext_bwa_index"] ? Channel.of (params.genomes[genome_name]["ext_bwa_index"]).first () : Channel.empty ()
		ch_ext_dict = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["ext_dict"] ? Channel.of (params.genomes[genome_name]["ext_dict"]).first () : Channel.empty ()
		ch_ext_fasta = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["ext_fasta"] ? Channel.of (params.genomes[genome_name]["ext_fasta"]).first () : Channel.empty ()
		ch_ext_fasta_index = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["ext_fasta_index"] ? Channel.of (params.genomes[genome_name]["ext_fasta_index"]).first () : Channel.empty ()
		ch_genome_base = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["genome_base"] ? Channel.of (params.genomes[genome_name]["genome_base"]).first () : Channel.empty ()
		ch_dir = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["dir"] ? Channel.of (params.genomes[genome_name]["dir"]).first () : Channel.empty ()
		ch_chrom_names = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["names"] && params.genomes[genome_name]["names"]["auto_sex"] ? Channel.value (params.genomes[genome_name]["names"]["auto_sex"]) : Channel.empty ()
		ch_chrom_names_auto = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["names"] && params.genomes[genome_name]["names"]["auto"] ? Channel.value (params.genomes[genome_name]["names"]["auto"]) : Channel.empty ()
		chrom_n = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["names"] && params.genomes[genome_name]["names"]["auto_sex"] ? params.genomes[genome_name]["names"]["auto_sex"].size () : 0

		cache_genome_url_bwa_index (genome_name, ch_genome_base, ch_ext_bwa_index)
		cache_genome_url_dict (genome_name, ch_genome_base, ch_ext_dict)
		cache_genome_url_fasta (genome_name, ch_genome_base, ch_ext_fasta)
		cache_genome_url_fasta_index (genome_name, ch_genome_base, ch_ext_fasta_index)

		ch_interval_list = ch_chrom_names.flatMap ().collectFile (name: 'interval_names.tsv', newLine: true)
		interval_bed (cache_genome_url_dict.out.result, ch_interval_list)

	emit:
		bwa_index        = cache_genome_url_bwa_index.out.result
		fasta            = cache_genome_url_fasta.out.result
		fasta_index      = cache_genome_url_fasta_index.out.result
		dict             = cache_genome_url_dict.out.result
		dir              = ch_dir
		chrom_names      = ch_chrom_names
		chrom_names_auto = ch_chrom_names_auto
		interval_bed     = interval_bed.out.result.first ()
		_chrom_n         = chrom_n
		
}

workflow GENOME_ANNOTATION
{
	take:
		genome_name

	main:
		if ( genome_name == null ) { exit 1, "[MoCaSeq] error: Genome name not found. Check params.genome_build." }
		ch_fasta_index_flat = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["fasta_index_flat"] ? Channel.of (params.genomes[genome_name]["fasta_index_flat"]).first () : Channel.empty ()
		ch_par_interval_bed = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["par_bed"] ? Channel.of (params.genome_annotations[genome_name]["par_bed"]).first () : Channel.empty ()
		ch_gc_wig = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["gc_wig"] ? Channel.of (params.genome_annotations[genome_name]["gc_wig"]) : Channel.empty ()
		ch_map_wig = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["map_wig"] ? Channel.of (params.genome_annotations[genome_name]["map_wig"]) : Channel.empty ()
		ch_common_vcf = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["common_vcf"] ? Channel.of (params.genome_annotations[genome_name]["common_vcf"]).first () : Channel.empty ()
		ch_gencode_genes_bed = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["gencode_genes_bed"] ? Channel.of (params.genome_annotations[genome_name]["gencode_genes_bed"]).first () : Channel.empty ()
		ch_micro_satellite = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["micro_satellite"] ? Channel.of (params.genome_annotations[genome_name]["micro_satellite"]).first () : Channel.empty ()

		cache_genome_url_common_vcf (genome_name, ch_common_vcf, Channel.value (["", "tbi"]))

		bash_expand_path_gc (ch_gc_wig)
		bash_expand_path_map (ch_map_wig)

	emit:
		fasta_index_flat = ch_fasta_index_flat
		par_interval_bed = ch_par_interval_bed
		gc_wig = bash_expand_path_gc.out.splitText ()
		map_wig = bash_expand_path_map.out.splitText ()
		common_vcf = cache_genome_url_common_vcf.out.result
		gencode_genes_bed = ch_gencode_genes_bed
		micro_satellite = ch_micro_satellite
}


