
workflow PREPARE_GENOME {

	take:
		genome_name

	main:
		ch_fasta = params.genomes[genome_name]["fasta"]
		ch_chrom_sizes = Channel.empty ()

	emit:
		fasta       = ch_fasta
		chrom_sizes = ch_chrom_sizes

}

