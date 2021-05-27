
process sam_to_fastq_paired {
	tag "${meta.sampleName}"

	input:
		val (reference)
		tuple val (meta), val (type), path (bam)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.R1.fastq.gz"), path ("${meta.sampleName}.${type}.R2.fastq.gz"), emit: result

	script:
	"""#!

java -Xmx${RAM}G -Dpicard.useLegacyParser=false -jar $picard_dir/picard.jar SamToFastq \\
	-INPUT $bam \\
	-FASTQ ${meta.sampleName}.${type}.R1.fastq.gz \\
	-SECOND_END_FASTQ ${meta.sampleName}.${type}.R2.fastq.gz \\
	-INCLUDE_NON_PF_READS true \\
	-VALIDATION_STRINGENCY LENIENT

	"""
}


