
params.gatk = [:]

process mutect_single {
	tag "${meta.sampleName}"

	input:
		val (genome_build)
		val (reference)
		tuple val (meta), val (type), path (bam)
		each (interval)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.m2.${interval}.vcf"), path ("${meta.sampleName}.${type}.m2.${interval}.f1r2.tar.gz"), emit: result

	script:
	"""#!/usr/bin/env bash

java -Xmx${params.gatk.ram}G -jar ${params.gatk.jar} Mutect2 \\
	--native-pair-hmm-threads 4 \\
	--reference ${reference} \\
	--intervals ${interval} \\
	--input ${bam} \\
	--f1r2-tar-gz ${meta.sampleName}.${type}.m2.${interval}.f1r2.tar.gz \\
	--output ${meta.sampleName}.${type}.m2.${interval}.vcf \\
	2> ${meta.sampleName}.${type}.m2.${interval}.log \\
	> ${meta.sampleName}.${type}.m2.${interval}.out

	"""

	stub:
	"""#!/usr/bin/env bash
#cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.m2.${interval}.vcf .
#cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.m2.${interval}.f1r2.tar.gz .
touch ${meta.sampleName}.${type}.m2.${interval}.vcf .
touch ${meta.sampleName}.${type}.m2.${interval}.f1r2.tar.gz .

	"""

}


process mutect_matched {
	 tag "${meta.sampleName}"

	input:
		val (genome_build)
		val (reference)
		tuple val (meta), path (bam_normal), path (bam_tumor)
		each (interval)

	output:
		tuple val (meta), val ("matched"), path ("${meta.sampleName}.matched.m2.${interval}.vcf"), path ("${meta.sampleName}.matched.m2.${interval}.f1r2.tar.gz"), emit: result

	script:
	"""#!/usr/bin/env bash

java -Xmx${params.gatk.ram}G -jar ${params.gatk.jar} Mutect2 \\
	--native-pair-hmm-threads 4 \\
	--reference ${reference} \\
	--intervals ${interval} \\
	--input ${bam_normal} \\
	--input ${bam_tumor} \\
	--normal-sample Normal --tumor-sample Tumor \\
	--f1r2-tar-gz ${meta.sampleName}.matched.m2.${interval}.f1r2.tar.gz \\
	--output ${meta.sampleName}.matched.m2.${interval}.vcf \\
	2> ${meta.sampleName}.matched.m2.${interval}.log \\
	> ${meta.sampleName}.matched.m2.${interval}.out

	"""

	stub:
	"""#!/usr/bin/env bash
#cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.matched.m2.${interval}.vcf .
#cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.matched.m2.${interval}.f1r2.tar.gz .
touch ${meta.sampleName}.matched.m2.${interval}.vcf .
touch ${meta.sampleName}.matched.m2.${interval}.f1r2.tar.gz .

	"""

}

process mutect_combine_vcf {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Mutect2", mode: "copy"

	input:
		val (genome_build)
		tuple val (meta), val (type), path ("*.interval.vcf"), path ("*.orientation_bias.tsv.gz")

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.m2.vcf.gz"), path ("${meta.sampleName}.${type}.m2.vcf.gz.tbi"), emit: result

	script:
	"""#!/usr/bin/env bash

vcf_files_first=\$(ls *.interval.vcf | head -n 1)
for f in *.interval.vcf;
do
	if [[ "\${f}" == "\${vcf_files_first}" ]]; then
		cat \${f}
	else
		cat \${f} | grep -v "^#"
	fi
done | bcftools sort --output-file ${meta.sampleName}.${type}.m2.vcf.gz --output-type z -
tabix -p vcf ${meta.sampleName}.${type}.m2.vcf.gz

cmd_learn_read_orientation="java -jar ${params.gatk.jar} LearnReadOrientationModel"
for f in *.orientation_bias.tsv.gz;
do
	cmd_learn_read_orientation="\${cmd_learn_read_orientation} --input \${f}"
done
cmd_learn_read_orientation="\${cmd_learn_read_orientation} --output ${meta.sampleName}.${type}.m2.read-orientation-model.tar.gz"

echo "\${cmd_learn_read_orientation}"
\${cmd_learn_read_orientation}

	"""

	stub:
	"""#!/usr/bin/env bash
#cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.m2.vcf.gz .
#cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/Mutect2/${meta.sampleName}.${type}.m2.vcf.gz.tbi .
touch ${meta.sampleName}.${type}.m2.vcf.gz
touch ${meta.sampleName}.${type}.m2.vcf.gz.tbi
	"""
}

