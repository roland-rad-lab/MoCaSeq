
params.strelka = [:]

process strelka_matched {

	tag "${meta.sampleName}"

	input:
		val (reference)
		val (interval_bed)
		tuple val (meta), path (bam_normal), path (bai_normal), path (bam_tumor), path (bai_tumor)
		path (candidate_small_indels_vcf)

	output:
		tuple val (meta), path("Strelka/results/variants/somatic.snvs.vcf.gz"), path ("Strelka/results/variants/somatic.indels.vcf.gz") emit: result

	script:
	"""#!/usr/bin/env bash
tabix -p bed ${interval_bed}
python2 ${params.strelka.dir}/bin/configureStrelkaSomaticWorkflow.py \\
	--normalBam ${bam_normal} \\
	--tumorBam ${bam_tumor} \\
	--ref ${reference} \\
	--runDir Strelka \\
	--indelCandidates ${candiate_small_indels_vcf} \\
	--callRegions ${interval_bed}

python2 Strelka/runWorkflow.py -m local -j ${params.strelka.threads}
	"""
}

