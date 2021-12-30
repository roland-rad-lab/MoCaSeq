
process mutect_filter_result_impact {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Mutect2", mode: "copy"

	input:
		val (genome_build)
		val (cgc_tsv)
		tuple val (meta), val (type), path (result_tsv)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.txt")
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.txt")
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt")

	script:
	"""#!/usr/bin/env Rscript

library (dplyr)

data <- read.table (file="${result_tsv}",sep="\\t",header=T,stringsAsFactors=F)
head (data)

data_cgc <- read.table (file="${cgc_tsv}",sep="\\t",header=T,stringsAsFactors=F)
head (data_cgc)

data_rare_impact <- data %>%
	dplyr::filter (ANN....IMPACT %in% c("HIGH", "MODERATE")) %>%
	data.frame

write.table (data_rare_impact %>% dplyr::mutate (dplyr::across (dplyr::everything (),~ ifelse (is.na (.x),"",.x))) %>% data.frame,file="${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.txt",sep="\\t",row.names=F,quote=F)

data_rare_impact_cgc <- data_rare_impact %>%
	dplyr::left_join (data_cgc %>% dplyr::select (Gene.Symbol) %>% data.frame,by=c("ANN....GENE"="Gene.Symbol")) %>%
	data.frame

write.table (data_rare_impact_cgc %>% dplyr::mutate (dplyr::across (dplyr::everything (),~ ifelse (is.na (.x),"",.x))) %>% data.frame,file="${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt",sep="\\t",row.names=F,quote=F)

	"""
}

process mutect_filter_result_impact_rare {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/Mutect2", mode: "copy"

	input:
		val (genome_build)
		val (cgc_tsv)
		val (tru_sight)
		tuple val (meta), val (type), path (result_tsv)

	output:
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.txt")
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.txt")
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt")
		tuple val (meta), val (type), path ("${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.CGC.TruSight.txt")

	script:
	"""#!/usr/bin/env Rscript

library (dplyr)

data <- read.table (file="${result_tsv}",sep="\\t",header=T,stringsAsFactors=F)
head (data)

data_cgc <- read.table (file="${cgc_tsv}",sep="\\t",header=T,stringsAsFactors=F)
head (data_cgc)

data_tru_sight <- read.table (file="${tru_sight}",sep="\\t",header=F,stringsAsFactors=F)
head (data_tru_sight)

data_rare <- data %>%
	dplyr::mutate (dplyr::across (AF, ~ dplyr::if_else (is.na (.x),0,.x))) %>%
	dplyr::mutate (dplyr::across (c(AC,AN), ~ dplyr::if_else (is.na (.x),0L,.x))) %>%
	dplyr::filter (G5=="false",AF < 0.1 & AN < 100 | AF <0.01 & AN >= 100) %>%
	data.frame

write.table (data_rare %>% dplyr::mutate (dplyr::across (dplyr::everything (),~ ifelse (is.na (.x),"",.x))) %>% data.frame,file="${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.txt",sep="\\t",row.names=F,quote=F)

data_rare_impact <- data_rare %>%
	dplyr::filter (ANN....IMPACT %in% c("HIGH", "MODERATE")) %>%
	data.frame

write.table (data_rare_impact %>% dplyr::mutate (dplyr::across (dplyr::everything (),~ ifelse (is.na (.x),"",.x))) %>% data.frame,file="${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.txt",sep="\\t",row.names=F,quote=F)

data_rare_impact_cgc <- data_rare_impact %>%
	dplyr::left_join (data_cgc %>% dplyr::select (Gene.Symbol) %>% data.frame,by=c("ANN....GENE"="Gene.Symbol")) %>%
	data.frame

write.table (data_rare_impact_cgc %>% dplyr::mutate (dplyr::across (dplyr::everything (),~ ifelse (is.na (.x),"",.x))) %>% data.frame,file="${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt",sep="\\t",row.names=F,quote=F)

data_rare_impact_cgc_tru_sight <- data_rare_impact_cgc %>%
	dplyr::left_join (data_tru_sight %>% dplyr::select (V1) %>% data.frame,by=c("ANN....GENE"="V1")) %>%
	data.frame

write.table (data_rare_impact_cgc_tru_sight %>% dplyr::mutate (dplyr::across (dplyr::everything (),~ ifelse (is.na (.x),"",.x))) %>% data.frame,file="${meta.sampleName}.${type}.Mutect2.NoCommonSNPs.OnlyImpact.CGC.TruSight.txt",sep="\\t",row.names=F,quote=F)

	"""

}

process mutect_signatures_matched {
	tag "${meta.sampleName}"

	input:
		val (genome_build)
		tuple val (meta), val (type), path (vcf), path (vcf_index)


	script:
	"""#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(SomaticSignatures))
suppressPackageStartupMessages(library(SomaticCancerAlterations))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(datasets))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(deconstructSigs))

BSgenome <- switch ("${meta.organism}",
	"human"={
		BSgenome.Hsapiens.UCSC.hg38
	},
	"mouse"={
		BSgenome.Mmusculus.UCSC.mm10
	},
	{
		stop ("Error: unrecognised organism '${meta.organism}'")
	})

system2 ("bcftools",args=c("query","--format","%CHROM\\\\t%POS\\\\t%REF\\\\t%ALT\\\\n","${vcf}"),stdout="${meta.sampleName}.${type}.Mutect2.tsv",wait=T)
data_sample <- read.table chr(file="${meta.sampleName}.${type}.Mutect2.tsv",sep="\\t",header=F,stringsAsFactors=F)
names (data_sample) <- c("chromosome", "pos", "ref", "alt")
head (data_sample)

sigs.input <- mut.to.sigs.input(mut.ref=data_sample, sample.id="${meta.sampleName}",chr="chromosome",pos="pos",ref="ref",alt="alt", bsg=BSgenome)

sigs.associated <- c(
	"Signature.1A",
	"Signature.2",
	"Signature.3",
	"Signature.4",
	"Signature.5",
	"Signature.6",
	"Signature.7",
	"Signature.8",
	"Signature.9",
	"Signature.10",
	"Signature.11",
	"Signature.12",
	"Signature.13",
	"Signature.14",
	"Signature.15",
	"Signature.16",
	"Signature.17",
	"Signature.18",
	"Signature.19",
	"Signature.20",
	"Signature.21"
)

sample <- whichSignatures (tumor.ref=sigs.input, signatures.ref=signatures.nature2013, associated=sigs.associataed, sample.id="${meta.sampleName}", contexts.needed=T, tri.counts.method="default", signature.cutoff=0.2)

pdf (file="${meta.sampleName}.${type}.Nature_Pie.pdf")
makePie (sample)
dev.off ()

pdf(file="${meta.sampleName}.${type}.Nature_Bar.pdf")
plotSignatures (sample)
dev.off ()

sample <- whichSignatures (tumor.ref=sigs.input, signatures.ref=signatures.cosmic, sample.id="${meta.sampleName}", contexts.needed=T, tri.counts.method="default", signature.cutoff=0.2)

pdf (file="${meta.sampleName}.${type}.Cosmic_Pie.pdf")
makePie (sample)
dev.off ()

pdf (file="${meta.sampleName}.${type}.Cosmic_Bar.pdf")
plotSignatures (sample)
dev.off ()


	"""

}




