
process dry_clean_detergent {

	publishDir "${params.output_base}/PON", mode: "copy", saveAs: { it.replaceFirst ("^PON/","${genome_build}.") }

	input:
		val (genome_build)
		val (intervals)
		val (par_region_bed)
		path ("*")
		path (normal_coverage_tsv)

	output:
		tuple path ("PON/normal_table.rds"), path ("PON/germline.markers.filtered.rds"), path ("PON/detergent.rds"), emit: result

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

grm <- readRDS("PON/germline.markers.rds")
grm_merged <- unlist(GenomicRanges::reduce (GenomicRanges::split (grm, ~germline.status)))
mcols(grm_merged)[,"germline.status.region"] <- names(grm_merged)

mcols(grm)[,"mi"] <- width(grm_merged)[GenomicRanges::findOverlaps (grm,grm_merged,select="first")]
mcols(grm)[mcols(grm)[,"mi"]==1000,"germline.status"] <- F

saveRDS(grm,file="PON/germline.markers.filtered.rds")

	"""
}

process dry_clean {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${meta.sampleName}/results/dryclean", mode: "copy"

	input:
		val (intervals)
		tuple path (normal_table_rds), path (germline_rds), path (detergent_rds)
		tuple val (meta), val (type), val(resolution), path (coverage_rds)

	output:
		tuple val (meta), val (type), val (resolution), path ("${meta.sampleName}.${type}.coverage.tsv"), emit: result
		tuple val (meta), val (type), val (resolution), path ("${meta.sampleName}.${type}.coverage.cnr"), emit: cnr

	script:
	"""#!/usr/bin/env Rscript
library (data.table)
library (IRanges)
library (GenomicRanges)
library (dryclean)
library (gUtils)

rescale_nn <- function (data)
{
	data_min = min (data,na.rm=T)
	data_max = max (data,na.rm=T)
	target_min=(data_max-data_min)/10000
	target_max=1

	return ( target_min + (data - data_min) * ((target_max - target_min)/(data_max-data_min)) )
}

intervals <- strsplit ("${intervals}", ",", fixed=T)[[1]]

coverage_data = readRDS("${coverage_rds}")
cov_out <- switch ("${type}",
	"Tumor"={
		start_wash_cycle(cov=coverage_data,detergent.pon.path="${detergent_rds}",whole_genome=T,chr=NA,germline.filter=T,germline.file="${germline_rds}",is.human=F,all.chr=intervals)	
	},
	"Normal"={
		start_wash_cycle(cov=coverage_data,detergent.pon.path="${detergent_rds}",whole_genome=T,chr=NA,germline.filter=T,germline.file="${germline_rds}",is.human=F,all.chr=intervals)
	},
	{
		stop ("Unsupported type: '${type}'")
	}
)

write.table (as.data.frame (cov_out),file="${meta.sampleName}.${type}.coverage.tsv",sep="\\t",quote=F,row.names=F)
cov_out_dt <- gUtils::gr2dt (cov_out)
cov_out_dt[,weight := abs(background.log) - median (background.log,na.rm=T)]
cov_out_dt[,weight := max(weight,na.rm=T) - weight]
cov_out_dt[,weight := rescale_nn (weight)]
cov_out_dt[is.na (log.reads),weight := cov_out_dt[,min(weight,na.rm=T)]]
cov_out_dt[,gene := ""]
setnames (cov_out_dt,c("seqnames","foreground.log","input.read.counts"),c("chromosome","log2","depth"))
write.table (cov_out_dt[,.(chromosome,start,end,gene,depth,log2,weight)],file="${meta.sampleName}.${type}.coverage.cnr",sep="\t",quote=F,row.names=F)

	"""
}

