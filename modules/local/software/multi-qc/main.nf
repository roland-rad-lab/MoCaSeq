


process qc {

	input:
		val (genome_build)
		path (analysis_dir)
		path (output_dir)


	script:
	"""
multiqc -o ${output_dir} --comment "MoCaSeq for ${genome_build} using analysis ${analysis_dir}" ${analysis_dir}
	"""
}

