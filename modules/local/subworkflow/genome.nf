
params.genomes = [:]

workflow PREPARE_GENOME
{

	take:
		genome_name

	main:
		ch_fasta = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["fasta"] ? Channel.fromPath (params.genomes[genome_name]["fasta"]) : Channel.empty ()
		ch_dict = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["dict"] ? Channel.fromPath (params.genomes[genome_name]["dict"]) : Channel.empty ()
		ch_chrom_names = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["names"] && params.genomes[genome_name]["names"]["auto_sex"] ? Channel.fromList (params.genomes[genome_name]["names"]["auto_sex"]) : Channel.empty ()
		chrom_n = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["names"] && params.genomes[genome_name]["names"]["auto_sex"] ? params.genomes[genome_name]["names"]["auto_sex"].size () : 0
		ch_chrom_sizes = Channel.empty ()

	emit:
		fasta       = ch_fasta
		dict        = ch_dict
		chrom_names = ch_chrom_names
		_chrom_n    = chrom_n
		chrom_sizes = ch_chrom_sizes
}

