
params.jabba = [:]

process jabba_matched {
	tag "${meta.sampleName}"


	script:
	"""#!/usr/bin/env Rscript
library(JaBbA)

system(paste("mkdir -p JaBbA"))

junctions =
coverage =

jab = JaBbA(
	## below are required positional arguments
	junctions = opt$junctions,
	coverage = opt$coverage,
	## below are junction related options
	juncs.uf = NULL,
	blacklist.junctions = NULL,
	whitelist.junctions = NULL,
	geno = F,
	indel = "exclude",
	cfield = opt$cfield,
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
	seg = "",
	max.na = 0.1,
	blacklist.coverage = NULL,
	nseg = "",
	hets = "",
	purity = scan(text = opt$purity, what = numeric(), sep = ",", quiet = T),
	ploidy = scan(text = opt$ploidy, what = numeric(), sep = ",", quiet = T),
	pp.method = "sequenza",
	## TODO: min.nbins, 5 by default
	cn.signif = 1e-5,
	## below are optimization related options
	slack = 100,
	loose.penalty.mode = "boolean",
	tilim = 6000,
	epgap = 0.01,
	## TODO use.gurobi = opt$gurobi,
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

}

