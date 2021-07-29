
process dry_clean_detergent {

	publishDir "${params.output_base}/PON", mode: "copy", saveAs: { it.replaceFirst ("^PON/","${genome_build}.") }

	input:
		val (genome_build)
		val (intervals)
		val (par_region_bed)
		path ("*")
		path (normal_coverage_tsv)

	output:
		tuple path ("PON/normal_table.rds"), path ("PON/germline.markers.rds"), path ("PON/detergent.rds"), emit: result

	script:
	"""#!/usr/bin/env Rscript

library (data.table)
library (IRanges)
library (GenomicRanges)
library (dryclean)

intervals <- strsplit ("${intervals}", ",", fixed=T)[[1]]

data_normal_coverage <- read.table (file="${normal_coverage_tsv}",sep="\\t",header=T,stringsAsFactors=F)
head (data_normal_coverage)
if ( nrow (data_normal_coverage) < 2 ) { stop ("You must supply more than one sample to calculate PON.") }

saveRDS (as.data.table (data_normal_coverage),file="normal_table.rds")

data_par <- read.table (file="${par_region_bed}",sep="\\t",header=F,stringsAsFactors=F)
names (data_par) <- c("seqnames", "start", "end")
head (data_par)
saveRDS (GenomicRanges::makeGRangesFromDataFrame (data_par),file="par.rds")

system("mkdir PON")


detergent <- prepare_detergent (normal.table.path="normal_table.rds",save.pon=T,path.to.save="PON",num.cores=1,use.all=T,PAR.file="par.rds",all.chr=intervals)
decomp_paths <- character (nrow (data_normal_coverage))

for ( i in seq_along (data_normal_coverage[,"sample"]) )
{
	sample_name <- data_normal_coverage[i,"sample"]
	sample_cov_path <- data_normal_coverage[i,"normal_cov"]

	decomp <- start_wash_cycle (cov=readRDS (sample_cov_path),detergent.pon.path="PON/detergent.rds",whole_genome=T,chr=NA,germline.filter=F,is.human=F,all.chr=intervals)
	head (decomp)

	sample_decomp_path <- paste ("PON/",sample_name,".decomp.rds",sep="")
	decomp_paths[i] <- sample_decomp_path
	saveRDS (decomp,file=sample_decomp_path)
}

saveRDS (as.data.table (cbind (data_normal_coverage,decomposed_cov=decomp_paths)),file="PON/normal_table.rds")

identify_germline (normal.table.path="PON/normal_table.rds",path.to.save="PON",save.grm=T,signal.thresh=0.5,pct.thresh=0.98,all.chr=intervals)

	"""
}

process dry_clean {
	tag "${meta.sampleName}"

	input:
		val (intervals)
		tuple path (normal_table_rds), path (germline_rds), path (detergent_rds)
		tuple val (meta), val (type), path (tumor_coverage_rds)

	output:
		tuple val (meta), path ("${meta.sampleName}.Tumor.coverage.tsv"), emit: result

	script:
	"""#!/usr/bin/env Rscript
library (data.table)
library (IRanges)
library (GenomicRanges)
library (dryclean)

intervals <- strsplit ("${intervals}", ",", fixed=T)[[1]]

coverage_data = readRDS("${tumor_coverage_rds}")
cov_out = start_wash_cycle(cov=coverage_data,detergent.pon.path="${detergent_rds}",whole_genome=T,chr=NA,germline.filter=T,germline.file="${germline_rds}",is.human=F,all.chr=intervals)

write.table (as.data.frame (cov_out),file="${meta.sampleName}.Tumor.coverage.tsv",sep="\\t",quote=F,row.names=F)

	"""
}

