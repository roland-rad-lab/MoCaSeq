
workflow PREPARE_GENOME {

	take:
		genome_name

	main:
		ch_fasta = params.genomes[genome_name]["fasta"]
		ch_chrom_sizes = Channel.empty ()

	emit:
		out.fasta       = ch_fasta
		out.chrom_sizes = ch_chrom_sizes

}

