
params.genomes = [:]

workflow PREPARE_GENOME
{

	take:
		genome_name

	main:
		ch_fasta = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["fasta"] ? Channel.fromPath (params.genomes[genome_name]["fasta"]) : Channel.empty ()
		ch_chrom_names = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["names"] && params.genomes[genome_name]["names"]["auto_sex"] ? Channel.fromList (params.genomes[genome_name]["names"]["auto_sex"]) : Channel.empty ()
		ch_chrom_sizes = Channel.empty ()

	emit:
		fasta       = ch_fasta
		chrom_names = ch_chrom_names
		chrom_sizes = ch_chrom_sizes
}

