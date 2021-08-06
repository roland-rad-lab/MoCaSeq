
include { qc } from "../software/multi-qc/main"

workflow COHORT_QC {

	take:
		genome_build
		ch_analysis_dir
		ch_output_dir

		ch_analysis_dir_verified = ch_analysis_dir.map {
			if ( !it.exists () ) exit 1, "[MoCaSeq] error: Analysis dir '${it}' does not exist."
			if ( !it.isDirectory () ) exit 1, "[MoCaSeq] error: Analysis dir '${it}' must be a directory."
			def analysis_dir_genome = it.resolve (genome_build)
			if ( !analysis_dir_genome.exists () ) println "[MoCaSeq] warning: Did not find '${analysis_dir_genome}'"
			analysis_dir_genome
		}
		ch_output_dir_verified = ch_output_dir.map {
			if ( !it.exists () ) exit 1, "[MoCaSeq] error: Output dir '${it}' does not exist."
			if ( !it.isDirectory () ) exit 1, "[MoCaSeq] error: Output dir '${it}' must be a directory."
			def output_dir_genome = it.resolve (genome_build)
			output_dir_genome.mkdir ()
			output_dir_genome
		}

	main:
		qc (genome_build, ch_analysis_dir_verified, ch_output_dir_verified)

}

