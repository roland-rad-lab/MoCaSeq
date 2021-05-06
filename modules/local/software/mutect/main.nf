
params.gatk = [:]

process mutect_matched {
	 tag "${meta.sampleName}"

	input:
		val (reference)
		tuple val (meta), path (bam_normal), path (bam_tumor)
		each (interval)

	output:
		tuple val (meta), path ("${meta.sampleName}.matched.m2.${interval}.vcf"), emit: vcf
		tuple val (meta), path ("${meta.sampleName}.matched.m2.${interval}.f1r2.tar.gz"), emit: orientation_bias

	script:
	"""#!/usr/bin/env bash

#java -Xmx${params.gatk.ram}G -jar ${params.gatk.jar} Mutect2 \\
#	--native-pair-hmm-threads 4 \\
#	--reference ${reference} \\
#	--intervals ${interval} \\
#	--input ${bam_normal} \\
#	--input ${bam_tumor} \\
#	--normal-sample Normal --tumor-sample Tumor \\
#	--f1r2-tar-gz ${meta.sampleName}.matched.m2.${interval}.f1r2.tar.gz \\
#	--output ${meta.sampleName}.matched.m2.${interval}.vcf \\
#	-bamout ${meta.sampleName}.matched.m2.${interval}.bam \\
#	--assembly-region-out ${meta.sampleName}.matched.m2.${interval}.assembly.txt \\
#	2> ${meta.sampleName}.matched.m2.${interval}.log \\
#	> ${meta.sampleName}.matched.m2.${interval}.out

touch ${meta.sampleName}.matched.m2.${interval}.vcf
touch ${meta.sampleName}.matched.m2.${interval}.f1r2.tar.gz
	"""
}

process mutect_combine_vcf {
	tag "${meta.sampleName}"

	input:
		tuple val (meta), path ("*.interval.vcf")

	output:
		tuple val (meta), path ("${meta.sampleName}.matched.m2.vcf.gz"), path ("${meta.sampleName}.matched.m2.vcf.gz.tbi"), emit: result

	script:
	"""#!/usr/bin/env bash
rm *.interval.vcf
cp /opt/fake_files/GX4VZC.Normal.m2.1.vcf 2.interval.vcf
cp /opt/fake_files/GX4VZC.Normal.m2.2.vcf 2.interval.vcf
cat *.interval.vcf | bcftools sort --output-file ${meta.sampleName}.matched.m2.${interval}.vcf.gz --output-type z -
tabix -p vcf ${meta.sampleName}.matched.m2.${interval}.vcf.gz
	"""
}

