
process igv_track_depth {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Tracks", mode: "copy"

	input:
		val (genome_build)
		val (intervals)
		tuple path (interval_bed), path (interval_bed_index)
		tuple val (meta), val (type), path (bam), path (bai)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.depth.raw.all.bigWig")

	script:
	"""#!/usr/bin/env bash

zcat ${interval_bed} | awk 'BEGIN{FS=OFS="\\t";}{print \$1,\$3-\$2;}' > genome.sizes

for chromosome in ${intervals};
do
	echo "fixedStep chrom=\${chromosome} start=1 step=1" >> ${meta.sampleName}.${type}.depth.raw.all.wig
	samtools depth -ar \${chromosome} ${bam} | cut -f 3 >> ${meta.sampleName}.${type}.depth.raw.all.wig
done
wigToBigWig ${meta.sampleName}.${type}.depth.raw.all.wig genome.sizes ${meta.sampleName}.${type}.depth.raw.all.bigWig
rm ${meta.sampleName}.${type}.depth.raw.all.wig

	"""
}

process igv_track_cnr {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Tracks", mode: "copy"

	input:
		val (genome_build)
		tuple path (interval_bed), path (interval_bed_index)
		tuple val (meta), val (type), val(coverage_source), path (cnr)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.${coverage_source}.bigWig")

	script:
	"""#!/usr/bin/env Rscript

library (dplyr)

data <- read.table (file="${cnr}",sep="\\t",header=T,stringsAsFactors=F)
data_genome_sizes <- read.table (file=gzfile ("${interval_bed}"),sep="\\t",header=F,stringsAsFactors=F)
names (data_genome_sizes) <- c("chromosome", "start", "end")

head (data)
head (data_genome_sizes)

write.table (data_genome_sizes %>% dplyr::select (chromosome,end) %>% data.frame ,file="genome.sizes",sep="\\t",row.names=F,col.names=F,quote=F)

data_bed <- switch ("${coverage_source}",
		{
			data
		}) %>%
		dplyr::filter (chromosome %in% !!data_genome_sizes[,"chromosome"]) %>%
		dplyr::select (chromosome,start,end,log2) %>%
		dplyr::mutate (start=as.numeric (start)) %>%
		dplyr::arrange (chromosome,start) %>%
		dplyr::mutate (start=format(start,scientific=F,trim=T)) %>%
		data.frame

write.table (data_bed,file="${meta.sampleName}.${type}.${coverage_source}.bedGraph",sep="\\t",row.names=F,col.names=F,quote=F)

system ("bedGraphToBigWig ${meta.sampleName}.${type}.${coverage_source}.bedGraph genome.sizes ${meta.sampleName}.${type}.${coverage_source}.bigWig")
file.remove ("${meta.sampleName}.${type}.${coverage_source}.bedGraph")

	"""
}

process igv_track_cns {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Tracks", mode: "copy"

	input:
		val (genome_build)
		tuple val (meta), val (type), val (coverage_source), path (cns)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.${coverage_source}.bedGraph")

	script:
	"""#!/usr/bin/env Rscript

library (dplyr)

data <- read.table (file="${cns}",sep="\\t",header=T,stringsAsFactors=F)
#head (data)

data_bed <- data %>%
	dplyr::select (chromosome,start,end,log2) %>%
	data.frame

write.table (data_bed,file="${meta.sampleName}.${type}.${coverage_source}.bedGraph",sep="\\t",quote=F,row.names=F)

	"""
}

process igv_track_rds {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Tracks", mode: "copy"

	input:
		val (genome_build)
		tuple path (interval_bed), path (interval_bed_index)
		val (coverage_source)
		tuple val (meta), val (type), val (resolution), path (coverage_rds)

	output:
		tuple val (meta), val (type), val (resolution), path ("${meta.sampleName}.${type}.${coverage_source}.${resolution}.bigWig")

	script:
	"""#!/usr/bin/env Rscript
library (dplyr)

interval_file <- gzfile ("${interval_bed}", 'rt')
data_interval <- read.table (file=interval_file,sep="\\t",header=F,stringsAsFactors=F)
names (data_interval) <- c("Chrom", "Start", "End")
head (data_interval)

write.table (data_interval %>% dplyr::mutate (Size=End-Start) %>% dplyr::select (Chrom,Size) %>% data.frame,file="intervals.sizes",sep="\\t",row.names=F,col.names=F,quote=F)

data <- as.data.frame (readRDS ("${coverage_rds}")) %>%
	dplyr::filter (!is.na(reads.corrected)) %>%
	dplyr::select (seqnames,start,end,reads.corrected) %>%
	dplyr::mutate (seqnames=as.character(seqnames)) %>%
	dplyr::arrange (seqnames,start) %>%
	data.frame

write.table (data,file="${meta.sampleName}.${type}.${coverage_source}.${resolution}.bedGraph",row.names=F,col.names=F,sep="\\t",quote=F)
system ("bedGraphToBigWig ${meta.sampleName}.${type}.${coverage_source}.${resolution}.bedGraph intervals.sizes ${meta.sampleName}.${type}.${coverage_source}.${resolution}.bigWig")
system ("rm ${meta.sampleName}.${type}.${coverage_source}.${resolution}.bedGraph")
	"""
}

process igv_track_vcf_sv_extract {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Tracks", mode: "copy"

	input:
		val (genome_build)
		tuple val (meta), val (type), path (vcf)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.vcf.tsv"), emit: result


	script:
	"""#!/usr/bin/env bash

echo "extract '${type}'"

case "${type}" in
	JaBbA)
	echo -e "SAMPLE\\tCHROM\\tPOS\\tID\\tREF\\tALT\\tINFO.MATEID\\tINFO.EVENT\\tINFO.SVTYPE\\tINFO.CNADJ\\tINFO.CNRADJ\\tINFO.CN\\tINFO.JUNCID\\tINFO.JABID\\tINFO.RJABID\\tGT\\tCN\\tRCN\\tSCN" | tr "[A-Z]" "[a-z]" > ${meta.sampleName}.${type}.vcf.tsv
	bcftools view -f PASS ${vcf} | bcftools query -f '[%SAMPLE\\t%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/MATEID\\t%INFO/EVENT\\t%INFO/SVTYPE\\t%INFO/CNADJ\\t%INFO/CNRADJ\\t%INFO/CN\\t%INFO/JUNCID\\t%INFO/JABID\\t%INFO/RJABID\\t%GT\\t%CN\\t%RCN\\t%SCN\\n]' | awk 'BEGIN{FS=OFS="\\t";}{\$1="${meta.sampleName}";print \$0;}' >> ${meta.sampleName}.${type}.vcf.tsv
	;;
	Manta)
	echo -e "SAMPLE\\tCHROM\\tPOS\\tID\\tREF\\tALT\\tINFO.IMPRECISE\\tINFO.SVTYPE\\tINFO.SVLEN\\tINFO.END\\tINFO.CIPOS\\tINFO.CIEND\\tINFO.CIGAR\\tINFO.MATEID\\tINFO.EVENT\\tINFO.HOMLEN\\tINFO.HOMSEQ\\tINFO.SVINSLEN\\tINFO.SVINSSEQ\\tINFO.LEFT_SVINSSEQ\\tINFO.RIGHT_SVINSSEQ\\tINFO.BND_DEPTH\\tINFO.MATE_BND_DEPTH\\tINFO.SOMATIC\\tINFO.SOMATICSCORE\\tINFO.JUNCTION_SOMATICSCORE\\tPR\\tSR" | tr "[A-Z]" "[a-z]" > ${meta.sampleName}.${type}.vcf.tsv
	bcftools view -f PASS ${vcf} | bcftools query -f '[%SAMPLE\\t%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/IMPRECISE\\t%INFO/SVTYPE\\t%INFO/SVLEN\\t%INFO/END\\t%INFO/CIPOS\\t%INFO/CIEND\\t%INFO/CIGAR\\t%INFO/MATEID\\t%INFO/EVENT\\t%INFO/HOMLEN\\t%INFO/HOMSEQ\\t%INFO/SVINSLEN\\t%INFO/SVINSSEQ\\t%INFO/LEFT_SVINSSEQ\\t%INFO/RIGHT_SVINSSEQ\\t%INFO/BND_DEPTH\\t%INFO/MATE_BND_DEPTH\\t%INFO/SOMATIC\\t%INFO/SOMATICSCORE\\t%INFO/JUNCTION_SOMATICSCORE\\t%PR\\t%SR\\n]' >> ${meta.sampleName}.${type}.vcf.tsv
	;;
	*)
	echo "No extraction implemented for '${type}'"
	exit 1
esac



	"""
}

process igv_track_vcf_sv {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Tracks", mode: "copy"

	input:
		val (genome_build)
		tuple val (meta), val (type), path (tsv)

	output:
		tuple val (meta), path ("${meta.sampleName}.${type}.bedpe")

	script:
	"""#!/usr/bin/env Rscript

library (purrr)
library (plyr)
library (dplyr)
library (stringr)

source ("${params.script_base}/VCF_Library.R")

data_sv_tsv <- read.table (file="${tsv}",sep="\\t",header=T,stringsAsFactors=F)
head (data_sv_tsv)

data_sv <- switch ("${type}",
	"JaBbA"={
		bedpe_from_jabba (data_sv_tsv)
	},
	{
		data
	}
) %>%
	dplyr::select (chrom1,start1,end1,chrom2,start2,end2,name,score,strand1,strand2) %>%
	data.frame

head (data_sv)
write.table (data_sv,file="${meta.sampleName}.${type}.bedpe",sep="\\t",row.names=F,col.names=F,quote=F)

	"""
}


