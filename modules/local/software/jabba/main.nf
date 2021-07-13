
params.jabba = [:]

process jabba_matched {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${meta.sampleName}/results/JaBbA", mode: "copy", saveAs: { it.replaceFirst ("^JaBbA/","") }

	input:
		tuple val (meta), path (manta_vcf), path (coverage_tsv), path (segments_tsv), val (purity), val (ploidy)

	output:
		tuple val (meta), path ("JaBbA/jabba.rds"), emit: result

	script:
	"""#!/usr/bin/env Rscript
options(error=function()traceback(2))

library(dplyr)
library(JaBbA)

system("rm -r JaBbA")
system("mkdir -p JaBbA")

data_coverage <- read.table (file="${coverage_tsv}",sep="\\t",header=T,stringsAsFactors=F)
data_segments <- read.table (file="${segments_tsv}",sep="\\t",header=T,stringsAsFactors=F)

junctions = "../${manta_vcf}"
coverage = data_coverage %>%
	mutate (log2Ratio=if_else(is.na (log2Ratio) | is.nan (log2Ratio),1,log2Ratio)) %>%
	dplyr::rename (seqnames=Chrom,start=Start,end=End) %>%
	dplyr::mutate (strand="*") %>%
	data.frame
segments = data_segments %>%
	dplyr::rename (seqnames=Chrom,start=Start,end=End) %>%
	dplyr::mutate (strand="*") %>%
	data.frame

head (segments)
purity = c(${purity})
ploidy = c(${ploidy})
cfield = "log2Ratio"

# Needed to edit JaBbA package
# In JaBbA/R/JaBbA.R:segstats
# utarget\$bad[is.nan(utarget\$mean)] = TRUE

jab = JaBbA(
	## below are required positional arguments
	junctions = junctions,
	coverage = coverage,
	## below are junction related options
	juncs.uf = NULL,
	blacklist.junctions = NULL,
	whitelist.junctions = NULL,
	geno = F,
	indel = "exclude",
	cfield = cfield,
	tfield = "tier",
	reiterate = 0,
	rescue.window = 1e3,
	nudge.balanced = T,
	## TODO: thresh.balanced, 500 default hardcoded
	edgenudge = 0.1,
	strict = F,
	all.in = F,
	## below are coverage related options
	field = "ratio",
	seg = segments,
	max.na = 0.1,
	blacklist.coverage = NULL,
	nseg = "",
	hets = "",
	purity = scan(text = purity, what = numeric(), sep = ",", quiet = T),
	ploidy = scan(text = ploidy, what = numeric(), sep = ",", quiet = T),
	pp.method = "sequenza",
	## TODO: min.nbins, 5 by default
	cn.signif = 1e-5,
	## below are optimization related options
	slack = 100,
	loose.penalty.mode = "boolean",
	tilim = 6000,
	epgap = 0.01,
	## TODO use.gurobi = opt\$gurobi,
	## TODO max.threads = Inf, max.mem = 16
	## below are general options
	outdir = "JaBbA",
	name = "${meta.sampleName}",
	mc.cores = ${params.jabba.threads},
	lp = F,
	ism = F,
	verbose = 2
	## TODO init, dyn.tuning
)

	"""

	stub:
	"""#!/usr/bin/env bash
mkdir -p JaBbA
cp ${params.stub_dir}/${meta.sampleName}/results/JaBbA/jabba.rds JaBbA/
	"""
}

process jabba_plot {
	tag "${meta.sampleName}"

	input:
		tuple val (meta), path (jabba_rds)

	script:
	"""#!/usr/bin/env Rscript

library (gGnome)

graph = gG (jabba="${jabba_rds}")
graph_events <- events (graph,verbose=F)
data_graph_events_type <- as.data.frame (graph_events\$meta\$event[, table (type)])

write.table (data_graph_events_type,file="${meta.sampleName}.events.counts.tsv",sep="\\t",quote=F,row.names=F)

for (i in seq_len (nrow (data_graph_events_type)) )
{
	event_type <- as.character (data_graph_events_type[i,"type"])
	event_count <- data_graph_events_type[i,"Freq"]
	plot_file_path <- paste ("${meta.sampleName}.event.",event_type,".pdf",sep="")

	pdf (file=plot_file_path,width=9)

	switch (event_type,
			dm={ e <- graph_events[dm>0];plot (e\$gt,streduce (e\$gr,1e6));title (paste(event_type," in ${meta.sampleName}")) },
			dup={ e <- graph_events[dup>0];plot (e\$gt,e\$footprint %>% GRanges %>% streduce(5e5));title (paste(event_type," in ${meta.sampleName}")) },
			inv={ e <- graph_events[simple>0];plot (e\$gt,e\$edges[grepl("^INV[0-9]+\$",simple)]\$shadow %>% streduce (1e5));title (paste(event_type," in ${meta.sampleName}")) },
			invdup={ e <- graph_events[simple>0];plot (e\$gt,e\$edges[grepl("^INVDUP[0-9]+\$",simple)]\$shadow %>% streduce (1e5));title (paste(event_type," in ${meta.sampleName}")) },
			tic={ e <- graph_events[tic>0];plot (e\$gt,e\$footprint %>% GRanges %>% streduce(5e5));title (paste(event_type," in ${meta.sampleName}")) },
			tra={ e <- graph_events[simple>0];plot (e\$gt,e\$edges[grepl("^TRA[0-9]+\$",simple)]\$shadow %>% streduce (1e5));title (paste(event_type," in ${meta.sampleName}")) },
			{
				stop (paste ("Event type '",event_type,"' is not implemented yet in sample ${meta.sampleName}",sep=""))
			}
	)

	dev.off ()
}

	"""

}


