
include { qc } from "../software/multi-qc/main"

workflow COHORT_QC {

	take:
		genome_build
		ch_output_dir

	main:

		def analysis_dir = file (params.output_base, glob: false)
		if ( !analysis_dir.exists () ) exit 1, "[MoCaSeq] error: Analysis dir '${analysis_dir}' does not exist."
		if ( !analysis_dir.isDirectory () ) exit 1, "[MoCaSeq] error: Analysis dir '${analysis_dir}' must be a directory."

		def analysis_genome_dir = analysis_dir.resolve (genome_build)
		if ( ! ( analysis_genome_dir.exists () && analysis_genome_dir.isDirectory () ) ) println "No analysis found for genome '${genome_build}'. The path '${analysis_genome_dir}' does not exist or was not a directory."

		ch_analysis_dir_verified = analysis_genome_dir.exists () && analysis_genome_dir.isDirectory () ? Channel.of (analysis_genome_dir) : Channel.empty ()

		ch_output_dir_verified = ch_output_dir.map {
			if ( !it.exists () ) exit 1, "[MoCaSeq] error: Output dir '${it}' does not exist."
			if ( !it.isDirectory () ) exit 1, "[MoCaSeq] error: Output dir '${it}' must be a directory."
			def output_dir_genome = it.resolve (genome_build)
			output_dir_genome.mkdir ()
			output_dir_genome
		}

		qc (genome_build, ch_analysis_dir_verified, ch_output_dir_verified)
}

