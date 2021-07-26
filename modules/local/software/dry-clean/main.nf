process dry_clean_detergent {

	publishDir "${params.output_base}/PON", mode: "copy"

	input:
		val (par_region_bed)
		path (normal_coverage_tsv)

	output:
		tuple path ("PON/normal_table.rds"), path ("PON/germline.markers.rds"), emit: result

	script:
	"""#!/usr/bin/env Rscript

library (data.table)
library (IRanges)
library (GenomicRanges)
library (dryclean)

data_normal_coverage <- read.table (file="${normal_coverage_tsv}",sep="\\t",header=T,stringsAsFactors=F)
head (data_normal_coverage)
if ( nrow (data_normal_coverage) < 2 ) { stop ("You must supply more than one sample to calculate PON.") }
dt_normal_coverage <- as.data.table (data_normal_coverage)
saveRDS (dt_normal_coverage,file="normal_table.rds")

data_par <- read.table (file="${par_region_bed}",sep="\\t",header=F,stringsAsFactors=F)
names (data_par) <- c("seqnames", "start", "end")
head (data_par)
saveRDS (GenomicRanges::makeGRangesFromDataFrame (data_par),file="par.rds")

system("mkdir PON")


detergent <- prepare_detergent (normal.table.path="normal_table.rds",save.pon=T,path.to.save="PON",num.cores=1,use.all=T,PAR.file="par.rds")

for ( i in seq_along (data_normal_coverage[,"sample"]) )
{
	sample_name <- data_normal_coverage[i,"sample"]
	sample_cov_path <- data_normal_coverage[i,"normal_cov"]

	decomp <- start_wash_cycle (cov=readRDS (sample_cov_path),detergent.pon.path="PON/detergent.rds",whole_genome=T,chr=NA,germline.filter=F)
	head (decomp)
}

grm = identify_germline (normal.table.path="normal_table.rds",path.to.save="PON",signal.thresh=0.5,pct.thresh=0.98)

	"""
}

process dry_clean {
	tag "${meta.sampleName}"

	input:
		tuple path (normal_table_rds), path (germline_rds)
		tuple val (meta), path (tumor_coverage_rds)

	script:
	"""#!/usr/bin/env

coverage_data = readRDS(tumor_coverage_rds)
cov_out = start_wash_cycle(cov=coverage_data,detergent.pon.path="PON",whole_genome=T,chr=NA,germline.filter=T,germline.file="germline.markers.rds")


write.table (data.frame (cov_out@unlistData),sep="\\t",quote=F,row.names=F)

	"""
}

