
params.gatk = [:]

process mutect_matched {
	 tag "${meta.sampleName}"

	input:
		val (reference)
		tuple val (meta), path (bam_normal), path (bam_tumor)
		each (interval)

	output:
		tuple val (meta), val ("matched"), path ("${meta.sampleName}.matched.m2.${interval}.vcf"), path ("${meta.sampleName}.matched.m2.${interval}.f1r2.tar.gz"), emit: result

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
		tuple val (meta), val (type), path ("*.interval.vcf"), path ("*.orientation_bias.tsv.gz")

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.m2.vcf.gz"), path ("${meta.sampleName}.${type}.m2.vcf.gz.tbi"), emit: result

	script:
	"""#!/usr/bin/env bash
rm *.interval.vcf
rm *.orientation_bias.tsv.gz
cp /opt/fake_files/GX4VZC.Normal.m2.1.vcf 2.interval.vcf
cp /opt/fake_files/GX4VZC.Normal.m2.2.vcf 1.interval.vcf
cp /opt/fake_files/GX4VZC.matched.1__186012250.1.m2.f1r2.tar.gz 2.orientation_bias.tsv.gz
cp /opt/fake_files/GX4VZC.matched.21__21344250.1.m2.f1r2.tar.gz 1.orientation_bias.tsv.gz

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
}

