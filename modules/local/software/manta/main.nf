
params.manta = [:]
params.snpeff = [:]
params.snpsift = [:]

process manta_matched {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Manta/variants", mode: "copy", pattern: "results/variants/*.vcf.gz*", saveAs: { it.replaceFirst ("^results/variants/","") }
	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Manta/evidence", mode: "copy", pattern: "results/evidence/*", saveAs: { it.replaceFirst ("^results/evidence/","") }
	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Manta/stats", mode: "copy", pattern: "results/stats/*", saveAs: { it.replaceFirst ("^results/stats/","") }

	input:
		val (genome_build)
		val (reference)
		tuple path (interval_bed), path (interval_bed_index)
		tuple val (meta), path (bam_normal), path (bai_normal), path (bam_tumor), path (bai_tumor)

	output:
		tuple val (meta), val ("matched"), path ("results/variants/somaticSV.vcf.gz"), path ("results/variants/somaticSV.vcf.gz.tbi"), emit: sv
		tuple val (meta), val ("matched"), path ("results/variants/candidateSmallIndels.vcf.gz"), path ("results/variants/candidateSmallIndels.vcf.gz.tbi"), emit: indel
		path ("results/evidence/*"), emit: evidence, optional: true
		path ("results/stats/*"), emit: stats, optional: true

	script:

	if ( "${meta.seqType}" == "wgs" )
	"""#!/usr/bin/env bash
python2 ${params.manta.dir}/bin/configManta.py \\
	--normalBam ${bam_normal} \\
	--tumorBam ${bam_tumor} \\
	--referenceFasta ${reference} \\
	--runDir . \\
	--callRegions ${interval_bed} \\
	--generateEvidenceBam

python2 runWorkflow.py -m local -j ${params.manta.threads}

	"""
	else if ( "${meta.seqType}" == "wex" )
	"""#!/usr/bin/env bash
python2 ${params.manta.dir}/bin/configManta.py \\
	--normalBam ${bam_normal} \\
	--tumorBam ${bam_tumor} \\
	--referenceFasta ${reference} \\
	--runDir . \\
	--callRegions ${interval_bed} \\
	--generateEvidenceBam \\
	--exome

python2 runWorkflow.py -m local -j ${params.manta.threads}

	"""
	else
		error "Invalid seqType: '${meta.seqType}' for sample: '${meta.sampleName}'"

	stub:
	"""#!/usr/bin/env bash
mkdir -p results/variants

if [[ "${params.stub_json_map?.manta_matched}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Manta/variants/somaticSV.vcf.gz results/variants/
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Manta/variants/somaticSV.vcf.gz.tbi results/variants/
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Manta/variants/candidateSmallIndels.vcf.gz results/variants/
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Manta/variants/candidateSmallIndels.vcf.gz.tbi results/variants/
fi

touch results/variants/somaticSV.vcf.gz
touch results/variants/somaticSV.vcf.gz.tbi
touch results/variants/candidateSmallIndels.vcf.gz
touch results/variants/candidateSmallIndels.vcf.gz.tbi

	"""
}

process manta_matched_post {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Manta", mode: "copy", pattern: "*.Manta.vcf.gz*"

	input:
		val (genome_build)
		val (snpeff_version)
		tuple val (meta), val(type), path (somatic_vcf), path (somatic_vcf_index)

	output:
		tuple val (meta), path ("${meta.sampleName}.Manta.annotated.vcf"), emit: result
		// It is only possible to publish process outputs
		tuple val (meta), path ("${meta.sampleName}.Manta.vcf.gz"), path ("${meta.sampleName}.Manta.vcf.gz.tbi"), emit: filtered_vcf

	script:

	"""#!/usr/bin/env bash

zcat ${somatic_vcf} | java -jar ${params.snpsift.jar} filter "( ( ( FILTER = 'PASS' ) & ( ( GEN[Tumor].SR[1] + GEN[Tumor].SR[0] ) * 0.05 <= GEN[Tumor].SR[1] ) ) | ( ( FILTER = 'PASS' ) & ( ( GEN[Tumor].PR[1] + GEN[Tumor].PR[0] ) * 0.05 <= GEN[Tumor].PR[1] ) ) )" | bgzip -c > ${meta.sampleName}.Manta.vcf.gz
tabix -p vcf ${meta.sampleName}.Manta.vcf.gz
bcftools stats ${meta.sampleName}.Manta.vcf.gz > ${meta.sampleName}.Manta.vcf.gz.stats

java -Xmx${params.snpeff.ram}g \\
	-jar ${params.snpeff.jar} \\
	${snpeff_version} -canon \\
	-csvStats ${meta.sampleName}.Manta.annotated.vcf.stats \\
	${meta.sampleName}.Manta.vcf.gz \\
	> ${meta.sampleName}.Manta.annotated.vcf

	"""

	stub:
	"""#!/usr/bin/env bash

if [[ "${params.stub_json_map?.manta_matched_post}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Manta/${meta.sampleName}.Manta.vcf.gz .
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Manta/${meta.sampleName}.Manta.vcf.gz.tbi .
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Manta/${meta.sampleName}.Manta.annotated.vcf .
fi

touch ${meta.sampleName}.Manta.vcf.gz
touch ${meta.sampleName}.Manta.vcf.gz.tbi
touch ${meta.sampleName}.Manta.annotated.vcf
	"""
}





