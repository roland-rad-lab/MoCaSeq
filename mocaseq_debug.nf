/* 
 * MoCaSeq debugging scripts
 * Here we execute snippets of the MoCaSeq pipeline for debugging
 * @author marcus.wagner@tum.de
 */


// write running log
log.info """\
    MoCaSeq nextflow test pipeline
    ===================================
    genomes_base 		: ${params.genomes_base}
    cache_base       	: ${params.cache_base}
    genome_build.human  : ${params.genome_build.human}
    """
    .stripIndent()


// helper for process cache_genome_url
boolean valid_uri (String maybe_url)
{
	try {
		java.net.URI.create (maybe_url);
		true
	} catch ( IllegalArgumentException | NullPointerException e ) {
		false
	}
}


// cache_genome_url process from modules/local/software/genome/main.nf
process cache_genome_url {
    input:
    val (genome_build)
    val (reference)
    val (extension_list)

    output:
    val (path_cached)

    exec:
	if (params.debug) { println "[MoCaSeq] debug: entered cache_genome_url process" }
	// if (params.debug) { println "[MoCaSeq] debug: genome_build:\n ${genome_build}\nreference:\n ${reference}\nextension_list:\n ${extension_list}" }
	def reference_local = reference;
	if ( extension_list.size () == 0 ) { exit 1, "[MoCaSeq] error: At least one file extension is required to append to '${reference_local}'" }
	if ( valid_uri (reference_local) )
	{
		def reference_uri = java.net.URI.create (reference_local)
		if ( reference_uri.scheme && reference_uri.scheme == "http" ) { exit 1, "[MoCaSeq] error: No support for http downloads use https, url requested was '${reference_local}'" }
		if ( reference_uri.scheme && reference_uri.scheme == "https" && reference_uri.path )
		{
			// println java.nio.file.Paths.get ("${params.cache_base}") // debug
			def cache_dir = java.nio.file.Paths.get ("${params.cache_base}").resolve ("genome_cache/${genome_build}");
			if (params.debug) { println "[MoCaSeq] debug: using cache_dir ${cache_dir}" }
			// println cache_dir // debug error in following line
			java.nio.file.Files.createDirectories (cache_dir);
			def file_name_base = reference_uri.path.tokenize ('/')[-1]
			path_cached = cache_dir.resolve ( [file_name_base, extension_list[0]].findAll { it }.join (".") ).toString ()

			extension_list.findAll { it != null }.each {
				def file_name_extended = [file_name_base, it].findAll { jt -> jt != "" }.join (".")
				if ( !valid_uri ( file_name_extended ) ) { exit 1, "[MoCaSeq] error: The extension '${file_name_extended}' was not a valid URI" }
				// download_uri (cache_dir.resolve (file_name_extended), reference_uri.resolve ( java.net.URI.create (file_name_extended)))
			}
		}
		else
		{
			path_cached = file (reference_local, glob: false).resolveSibling ( [file (reference_local, glob: false).fileName, extension_list[0]].findAll { it }.join (".") ).toString ()
		}
	}
	else
	{
		path_cached = file (reference_local, glob: false).resolveSibling ( [file (reference_local, glob: false).fileName, extension_list[0]].findAll { it }.join (".") ).toString ()
	}
}

// model workflow from modules/local/subworkflow/genome.nf:46,50,56
workflow PREPARE_GENOME {
	take:
		genome_name
	main:
		ch_ext_bwa_index = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["ext_bwa_index"] ? Channel.of (params.genomes[genome_name]["ext_bwa_index"]).first () : Channel.empty ()
		ch_genome_base = params.genomes && params.genomes[genome_name] && params.genomes[genome_name]["genome_base"] ? Channel.of (params.genomes[genome_name]["genome_base"]).first () : Channel.empty ()
		cache_genome_url (genome_name, ch_genome_base, ch_ext_bwa_index)
}

// workflow from main.nf
workflow HUMAN_MAP {
	if (params.debug) { println "[MoCaSeq] debug: entered HUMAN_MAP worfklow" }
	PREPARE_GENOME (params.genome_build.human)
}

