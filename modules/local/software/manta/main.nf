
params.manta = [:]

process interval_bed {

	input:
		val (dict)
		path (interval_list)

	output:
		path ("intervals.bed.gz"), emit: result

	script:
	"""#!/usr/bin/env Rscript

library (dplyr)
library (stringr)

dict_file_path = "${dict}"
intervals_file_path = "interval_names.tsv"
output_bed_file_path = "intervals.bed.gz"

data_dict <- read.table (file=dict_file_path,sep="\\t",stringsAsFactors=F,header=F,skip=1)
names (data_dict) <- c("line_type", "sequence_name_raw", "sequence_length_raw", "md5", "url")
head (data_dict)

data_seq_lengths <- data_dict %>%
  mutate (sequence_name=stringr::str_split_fixed (sequence_name_raw,":",2)[,2]) %>%
  mutate (sequence_length=as.numeric (stringr::str_split_fixed (sequence_length_raw,":",2)[,2])) %>%
  select (sequence_name,sequence_length) %>%
  data.frame

head (data_seq_lengths)

data_intervals <- read.table (file=intervals_file_path,sep="\\t",stringsAsFactors=F,header=F)
names (data_intervals) <- c("sequence_name")

data_output <- data_intervals %>%
	inner_join (data_seq_lengths,by="sequence_name") %>%
	mutate (start=1) %>%
	select (sequence_name,start,sequence_length) %>%
	data.frame

output_bed_file <-pipe (paste ("bgzip -c >",output_bed_file_path), "w")
write.table (data_output,file=output_bed_file,sep="\\t",row.names=F,col.names=F,quote=F)
close (output_bed_file)
"""

}

process manta_matched {
	tag "${meta.sampleName}"

	input:
		val (reference)
		val (interval_bed)
		tuple val (meta), path (bam_normal), path (bai_normal), path (bam_tumor), path (bai_tumor)

	output:
		tuple val (meta), path ("results/variants/somaticSV.vcf.gz"), emit: result

	script:
	"""#!/usr/bin/env bash
tabix -p bed ${interval_bed}
python2 ${params.manta.dir}/bin/configManta.py \\
	--normalBam ${bam_normal} \\
	--tumorBam ${bam_tumor} \\
	--referenceFasta ${reference} \\
	--runDir . \\
	--callRegions ${interval_bed} \\
	--generateEvidenceBam

python2 runWorkflow.py -m local -j ${params.manta.threads}

	"""
}

process manta_matched_post {
	tag "${meta.sampleName}"

	input:
		tuple val (meta), path (somatic_vcf)


	script:
	"""#!/usr/bin/env bash

zcat ${somatic_vcf} | java -jar ${params.annotation.snpeff.dir}/SnpSift.jar filter "( ( ( FILTER = 'PASS' ) & ( ( GEN[Tumor].SR[1] + GEN[Tumor].SR[0] ) * 0.05 <= GEN[Tumor].SR[1] ) ) | ( ( FILTER = 'PASS' ) & ( ( GEN[Tumor].PR[1] + GEN[Tumor].PR[0] ) * 0.05 <= GEN[Tumor].PR[1] ) ) )" | bgzip -c > ${meta.sampleName}.Manta.vcf.gz
tabix -p vcf ${meta.sampleName}.Manta.vcf.gz
bcftools stats ${meta.sampleName}.Manta.vcf.gz > ${meta.sampleName}.Manta.vcf.gz.stats

java -Xmx${params.annotation.snpeff.ram}g \\
	-jar ${params.annotation.snpeff.dir}/snpEff.jar \\
	${params.annotation.snpeff[${meta.organism}]} -canon \\
	-csvStats ${meta.sampleName}.Manta.annotated.vcf.stats \\
	${meta.sampleName}.Manta.vcf.gz \\
	> ${meta.sampleName}.Manta.annotated.vcf

	"""
}





