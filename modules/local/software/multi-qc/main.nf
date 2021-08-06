


process qc {

	input:
		val (genome_build)
		path (analysis_dir, stageAs: "analysis_dir")
		path (output_dir, stageAs: "output_dir")


	script:
	"""
multiqc -o ${output_dir} --comment "MoCaSeq for ${genome_build}" ${analysis_dir}
	"""
}

