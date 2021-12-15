#!/usr/bin/env nextflow

include {
	interval_bed
	cache_genome_url as cache_genome_url_bwa_index
	cache_genome_url as cache_genome_url_all_vcf
	cache_genome_url as cache_genome_url_common_vcf
	cache_genome_url as cache_genome_url_dict
	cache_genome_url as cache_genome_url_fasta
	cache_genome_url as cache_genome_url_fasta_index
	cache_genome_url as cache_genome_url_gc_wig
	cache_genome_url as cache_genome_url_gencode_genes_bed
	cache_genome_url as cache_genome_url_map_wig
	cache_genome_url as cache_genome_url_micro_satellite
	cache_genome_url as cache_genome_url_ref_flat
	cache_genome_url as cache_genome_url_sift_dbsnp
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
		ch_all_vcf = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["all_vcf"] ? Channel.of (params.genome_annotations[genome_name]["all_vcf"]).first () : Channel.empty ()
		ch_ref_flat = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["ref_flat"] ? Channel.of (params.genome_annotations[genome_name]["ref_flat"]).first () : Channel.empty ()
		ch_par_interval_bed = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["par_bed"] ? Channel.of (params.genome_annotations[genome_name]["par_bed"]).first () : Channel.empty ()
		ch_gc_wig = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["gc_wig"] ? Channel.of (params.genome_annotations[genome_name]["gc_wig"]) : Channel.empty ()
		ch_map_wig = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["map_wig"] ? Channel.of (params.genome_annotations[genome_name]["map_wig"]) : Channel.empty ()
		ch_common_vcf = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["common_vcf"] ? Channel.of (params.genome_annotations[genome_name]["common_vcf"]).first () : Channel.empty ()
		ch_gencode_genes_bed = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["gencode_genes_bed"] ? Channel.of (params.genome_annotations[genome_name]["gencode_genes_bed"]).first () : Channel.empty ()
		ch_micro_satellite = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["micro_satellite"] ? Channel.of (params.genome_annotations[genome_name]["micro_satellite"]).first () : Channel.empty ()
		ch_sift_dbsnp = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["sift_dbsnp"] ? Channel.of (params.genome_annotations[genome_name]["sift_dbsnp"]).first () : Channel.empty ()
		ch_sift_fields = params.genome_annotations && params.genome_annotations[genome_name] && params.genome_annotations[genome_name]["sift_fields"] ? Channel.of (params.genome_annotations[genome_name]["sift_fields"]).first () : Channel.empty ()

		ch_gc_wig_branched = ch_gc_wig.branch {
			uri: it.startsWith ("https://")
			other: true
		}
		ch_map_wig_branched = ch_map_wig.branch {
			uri: it.startsWith ("https://")
			other: true
		}

		cache_genome_url_all_vcf (genome_name, ch_all_vcf, Channel.value (["", "tbi"]))
		cache_genome_url_common_vcf (genome_name, ch_common_vcf, Channel.value (["", "tbi"]))
		cache_genome_url_gc_wig (genome_name, ch_gc_wig_branched.uri, Channel.value ([""]))
		cache_genome_url_gencode_genes_bed (genome_name, ch_gencode_genes_bed, Channel.value ([""]))
		cache_genome_url_map_wig (genome_name, ch_map_wig_branched.uri, Channel.value ([""]))
		cache_genome_url_micro_satellite (genome_name, ch_micro_satellite, Channel.value ([""]))
		cache_genome_url_ref_flat (genome_name, ch_ref_flat, Channel.value ([""]))
		cache_genome_url_sift_dbsnp (genome_name, ch_sift_dbsnp, Channel.value ([""]))

		bash_expand_path_gc (ch_gc_wig_branched.other)
		bash_expand_path_map (ch_map_wig_branched.other)

		ch_sift_sources = cache_genome_url_sift_dbsnp.out.result

	emit:
		par_interval_bed = ch_par_interval_bed
		all_vcf = cache_genome_url_all_vcf.out.result
		common_vcf = cache_genome_url_common_vcf.out.result
		gc_wig = bash_expand_path_gc.out.splitText ().mix (cache_genome_url_gc_wig.out.result)
		gencode_genes_bed = cache_genome_url_gencode_genes_bed.out.result
		map_wig = bash_expand_path_map.out.splitText ().mix (cache_genome_url_map_wig.out.result)
		micro_satellite = cache_genome_url_micro_satellite.out.result
		ref_flat = cache_genome_url_ref_flat.out.result
		sift_fields = ch_sift_fields
		sift_sources = ch_sift_sources
}


