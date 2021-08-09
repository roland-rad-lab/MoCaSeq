
process delly_matched_call {
	tag "${meta.sampleName}"

	input:
		val (genome_build)
		val (reference)
		tuple val(meta), path (bam_normal), path (bai_normal), path (bam_tumor), path (bai_tumor)

	output:
		tuple val(meta), path("${meta.sampleName}.pre.bcf"), path ("${meta.sampleName}.pre.bcf.csi"), emit: result

	script:
	"""#!/usr/bin/env bash
delly call \\
	-o ${meta.sampleName}.pre.bcf \\
	-g ${reference} \\
	${bam_tumor} \\
	${bam_normal}

	"""

	stub:
	"""#!/usr/bin/env bash
#cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Delly/${meta.sampleName}.pre.bcf .
#cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Delly/${meta.sampleName}.pre.bcf.csi .
touch ${meta.sampleName}.pre.bcf
touch ${meta.sampleName}.pre.bcf.csi
	"""
}

process delly_matched_filter {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Delly", mode: "copy"

	input:
		val (genome_build)
		tuple val(meta), path (delly_pre_bcf), path (delly_pre_bcf_index)

	output:
		tuple val(meta), path("${meta.sampleName}.delly.vcf.gz"), path ("${meta.sampleName}.delly.vcf.gz.tbi"), emit: result

	script:
	"""#!/usr/bin/env bash
cat <<"EOF" > Samples.tsv
Tumor	tumor
Normal	control
EOF

delly filter \\
	-f somatic \\
	-o ${meta.sampleName}.delly.bcf \\
	-s Samples.tsv \\
	${delly_pre_bcf}

bcftools view ${meta.sampleName}.delly.bcf -O z -o ${meta.sampleName}.delly.vcf.gz
tabix -p vcf ${meta.sampleName}.delly.vcf.gz
	"""

	stub:
	"""#!/usr/bin/env bash
#cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Delly/${meta.sampleName}.delly.vcf.gz .
#cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Delly/${meta.sampleName}.delly.vcf.gz.tbi .
touch ${meta.sampleName}.delly.vcf.gz
touch ${meta.sampleName}.delly.vcf.gz.tbi
	"""

}

