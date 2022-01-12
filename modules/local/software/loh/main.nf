
process loh_matched {
	tag "${meta.sampleName}"

	input:
		val (genome_build)
		tuple path (interval_bed), path (interval_bed_index)
		tuple val (meta), path (normal_vcf), path (normal_vcf_index), path (tumor_vcf), path (tumor_vcf_index)

	output:
		tuple val (meta), path ("${meta.sampleName}.VariantsForLOHGermline.tsv.gz"), emit: result

	script:
	"""#!/usr/bin/env bash
source ${params.script_base}/file_handling.sh
temp_file_b=\$(moc_mktemp_file . bed)
trap "rm \${temp_file_b}" EXIT

extract_if_zip ${interval_bed} interval_bed_extracted \${temp_file_b}

bcftools merge \\
	--regions-file \${interval_bed_extracted} \\
	--merge all \\
	--output-type u \\
	${normal_vcf} \\
	${tumor_vcf} \\
	| bcftools view \\
	-m2 -M2 \\
	--types snps,indels \\
	--output-type u \\
	| bcftools query \\
	--include "F_PASS(INFO/MMQ[1]>=60 & FORMAT/DP>=10) == 1.0" \\
	--format '[%SAMPLE\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%AF\\t%AD{0}\\t%AD{1}\\n]' \\
	| gzip > ${meta.sampleName}.VariantsForLOHGermline.tsv.gz
	"""

	stub:
	"""#!/usr/bin/env bash
touch ${meta.sampleName}.VariantsForLOHGermline.tsv.gz
	"""

}

process loh_matched_assign_alleles {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/LOH", mode: "copy"

	input:
		val (genome_build)
		val (reference)
		val (reference_index)
		val (intervals)
		tuple val (meta), path (loh_variants_tsv)

	output:
		tuple val (meta), path ("${meta.sampleName}.VariantsForLOH.txt"), emit: result

	script:
	"""#!/usr/bin/env python3.7
import gzip
import itertools
import collections

import pysam

reference_fasta = pysam.FastaFile ("${reference}")
chrom_order = collections.defaultdict (int, {v:k for k,v in enumerate ("${intervals}".split (","))})

top_fwd = set ([ 
	("A","C"),
	("A","G")
	])
top_rev = set ([
	("C","A"),
	("G","A")
	])
bot_fwd = set ([
	("T","C"),
	("T","G")
	])
bot_rev = set ([
	("C","T"),
	("G","T")
	])

unambiguous_top = top_fwd | top_rev
unambiguous_bot = bot_fwd | bot_rev

ambiguous_top = top_fwd | bot_fwd
ambiguous_bot = top_rev | bot_rev

header = "sample,chrom,pos,ref,alt,af,ad_ref,ad_alt".split (",")
output_header = "UniquePos,Chrom,Pos,Ref,Alt,Normal_Freq,Normal_RefCount,Normal_AltCount,Tumor_Freq,Tumor_RefCount,Tumor_AltCount,Genotype,strand_vector,Plot_Freq,genotype_ab,resolution,context".split (",")

def variant_key (data):
	return ("_".join ([data[k] for k in ["chrom", "pos", "ref", "alt"]]))

def position_key (data):
	return (tuple ([chrom_order[data[1]], int (data[2])]))

with open ("${meta.sampleName}.VariantsForLOH.txt", "w") as output_file:
	output = []
	output_file.write ("\\t".join (output_header))
	output_file.write ("\\n")
	with gzip.open ("${loh_variants_tsv}", "rt") as input_file:
		data = []
		for line in input_file:
			ldata = line.rstrip ().split ("\\t")
			data.append ({k:v for k,v in zip (header,ldata)})
		for k,g in itertools.groupby (sorted(data,key=variant_key),key=variant_key):
			lg = list (g)
			#print (k)
			sample_normal = next (filter (lambda x: x["sample"] == "Normal",lg),None)
			sample_tumor = 	next (filter (lambda x: x["sample"] == "Tumor",lg),None)
			if sample_normal == None or sample_tumor == None:
				raise Exception ("Failed to find both Normal and Tumor sample for %s" % k)
			#print (sample_normal)
			#print (sample_tumor)

			if float (sample_normal["af"]) > 0.7 or float (sample_normal["af"]) < 0.3:
				continue
			genotype = (sample_normal["ref"], sample_normal["alt"])
			genotype_ab = (sample_normal["ref"], sample_normal["alt"])
			strand = "NA"
			plot_freq = 0.0

			af_tumor = float (sample_tumor["af"])
			resolution = ""
			dbg_left = []
			dbg_right = []

			if len (genotype[0]) > 1 or len (genotype[1]) > 1:
				plot_freq = af_tumor
			elif genotype in unambiguous_top and genotype[0] == "A":
				plot_freq = af_tumor
				strand = "TOPua"
			elif genotype in unambiguous_top:
				plot_freq = 1.0 - af_tumor
				strand = "TOPua"
				genotype_ab = tuple (reversed(genotype))
			elif genotype in unambiguous_bot and genotype[0] == "T":
				plot_freq = 1.0 - af_tumor
				strand = "BOTTOMua"
			elif genotype in unambiguous_bot:
				plot_freq = af_tumor
				strand = "BOTTOMua"
				genotype_ab = tuple (reversed(genotype))
			else:
				# Resolve an ambiguous variant
				#print ("resolving")
				sequence_context = reference_fasta.fetch (sample_normal["chrom"],start=int (sample_normal["pos"])-101,end=int(sample_normal["pos"])+100)
				#print (sequence_context)
				#print (len(sequence_context))
				#print (sequence_context[100])
				#print (sequence_context[98:103])
				for scd in map (tuple,zip(sequence_context[:100][::-1],sequence_context[101:])):
					#print (scd)
					dbg_left.append (scd[0])
					dbg_right.append (scd[1])
					resolution = "".join (scd)
					if scd in ambiguous_top:
						strand = "TOP"
						plot_freq = af_tumor
						break
					if scd in ambiguous_bot:
						strand = "BOTTOM"
						plot_freq = 1.0 - af_tumor
						genotype_ab = tuple (reversed(genotype))
						break
			odata = []
			odata.append (k)
			odata.append (sample_normal["chrom"])
			odata.append (sample_normal["pos"])
			odata.append (sample_normal["ref"])
			odata.append (sample_normal["alt"])
			odata.append (sample_normal["af"])
			odata.append (sample_normal["ad_ref"])
			odata.append (sample_normal["ad_alt"])
			odata.append (sample_tumor["af"])
			odata.append (sample_tumor["ad_ref"])
			odata.append (sample_tumor["ad_alt"])
			odata.append ("".join (genotype))
			odata.append (strand)
			odata.append ("%0.3f" % plot_freq)
			odata.append ("".join (genotype_ab))
			odata.append (resolution)
			odata.append ("".join (list(reversed(dbg_left)) + ["@"] + dbg_right))
			output.append (odata)

	for odata in sorted (output, key=position_key):
		output_file.write ("\\t".join (odata))
		output_file.write ("\\n")

	"""

	stub:
	"""#!/usr/bin/env bash
if [[ "${params.stub_json_map?.loh_matched_assign_alleles}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/LOH/${meta.sampleName}.VariantsForLOH.txt .
fi

touch ${meta.sampleName}.VariantsForLOH.txt
	"""

}

process loh_matched_plot {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/LOH", mode: "copy"

	input:
		val (genome_build)
		tuple path (interval_bed), path (interval_bed_index)
		tuple val (meta), path (loh_tsv)

	output:
		path ("${meta.sampleName}.*.pdf")

	script:
	"""#!/usr/bin/env Rscript
library (dplyr)
library (ggplot2)
library (gridExtra)

interval_file <- gzfile ("${interval_bed}", 'rt')
data_interval <- read.table (file=interval_file,sep="\\t",header=F,stringsAsFactors=F)
names (data_interval) <- c("Chrom", "Start", "End")
head (data_interval)

data <- read.table (file="${loh_tsv}",sep="\\t",header=T,stringsAsFactors=F)
head (data)

data_interval_plot <- data_interval %>%
	dplyr::mutate (End=as.numeric (End)) %>%
	dplyr::mutate (CumulativeStart=cumsum (End)-End) %>%
	dplyr::mutate (CumulativeMidpoint=CumulativeStart+((Start+End)/2)) %>%
	data.frame

data_plot <- data %>%
	dplyr::inner_join (data_interval_plot,by="Chrom",suffix=c("",".Chrom")) %>%
	dplyr::mutate (Start.Genome=Pos+CumulativeStart) %>%
	data.frame

plot_types <- setNames (c("Plot_Freq","Tumor_Freq","Normal_Freq"),c("adjusted","raw","germline"))

for ( i in seq_along (plot_types) )
{
	g <- ggplot (data_plot) +
		geom_point (aes_string(x="Start.Genome",y=plot_types[i]),shape=".",colour="#ff80c3") +
		geom_vline (data=data_interval_plot,aes(xintercept=CumulativeStart)) +
		geom_text (data=data_interval_plot,aes(x=CumulativeMidpoint,y=1.1,label=Chrom),size=2) +
		coord_cartesian (ylim=c(0,1),expand=F,clip="off") +
		theme_bw () +
		theme (
			panel.grid.major.x=element_blank (),
			panel.grid.minor.x=element_blank (),
			axis.ticks.x=element_blank (),
			axis.text.x=element_blank (),
			plot.margin = unit(c(1,0.5,0.5,0.5), "cm")
		)
	pdf (file=paste ("${meta.sampleName}.LOH.",names (plot_types)[i],".genome.pdf",sep=""),width=16,height=4)
	print (g)
	dev.off ()

	chromosomes <- data_interval %>% pull (Chrom)
	plot_list <- vector ("list",length (chromosomes))

	for ( j in seq_along (chromosomes) )
	{
		plot_list[[j]] <- ggplot (data_plot %>% filter (Chrom==!!chromosomes[[j]]) %>% data.frame) +
			geom_point (aes_string(x="Pos",y=plot_types[i]),shape=".",colour="#ff80c3") +
			ylim (0,1) +
			labs (title=chromosomes[[i]]) +
			theme_bw () +
			theme (
				panel.grid.major.x=element_blank (),
				panel.grid.minor.x=element_blank (),
				axis.text.x=element_blank ()
			)
	}

	# Need to do this outside of pdf call to prevent blank first page
	p <- marrangeGrob (plot_list,nrow=1,ncol=1)

	pdf (file=paste ("${meta.sampleName}.LOH.",names (plot_types)[i],".chromosomes.pdf",sep=""),width=9)
	print (p)
	dev.off ()

}

	"""

	stub:
	"""#!/usr/bin/env bash
touch ${meta.sampleName}.LOH.adjusted.chromosomes.pdf
touch ${meta.sampleName}.LOH.adjusted.genome.pdf
	"""

}

process loh_matched_segment {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/LOH", mode: "copy"

	input:
		val (genome_build)
		val (intervals)
		tuple path (interval_bed), path (interval_bed_index)
		val (alt_haplotype)
		val (centromere)
		val (mappability)
		tuple val (meta), path (loh_tsv)

	output:
		tuple val (meta), path ("${meta.sampleName}.LOH.snp_filter.tsv"), path ("${meta.sampleName}.LOH.segs.tsv"), emit: result

	script:
	"""#!/usr/bin/env Rscript

library (dplyr)
library (GenomicRanges)
library (PSCBS)

# Following method by NikdAK

getLOHSegments <- function (data_gr)
{
	data_points <- as.data.frame (data_gr) %>%
		dplyr::rename (chromosome=seqnames,x=start,y=Plot_Freq) %>%
		dplyr::select (chromosome,x,y) %>%
		data.frame

	#print ("data_points")
	#print (head (data_points))

	data_points_dup <- data_points %>%
		dplyr::mutate (y=abs(y)) %>%
		dplyr::mutate (chromosome=as.character (chromosome)) %>%
		data.frame

	#print (head (data_points_dup))

	gaps <- PSCBS::findLargeGaps (data_points_dup, minLength = 2e+06) # smallest centromere is 2100000
	knownSegments <- PSCBS::gapsToSegments (gaps)
	fit <- PSCBS::segmentByCBS (data_points_dup, knownSegments = knownSegments, verbose = -10)

	#print ("fit")
	#print (fit)

	fit_gr <- as.data.frame (fit) %>%
		dplyr::filter (!is.na(mean)) %>%
		dplyr::mutate (mid=(start+end)/2) %>%
		makeGRangesFromDataFrame (keep.extra.columns=T)

	#print (fit_gr)

	fit_gr_data_gr.hits <- GenomicRanges::findOverlaps (fit_gr,data_gr)
	fit_gr_data_gr.hits_data <- cbind (data.frame (fit_gr_data_gr.hits), data.frame (LOHscore=mcols (data_gr)[subjectHits (fit_gr_data_gr.hits),"Plot_Freq"])) %>%
		dplyr::group_by (queryHits) %>%
		dplyr::summarise (mode=abs(density (LOHscore)\$x[which.max(density(LOHscore)\$y)])) %>%
		data.frame

	#print (head (fit_gr_data_gr.hits_data))
	mcols (fit_gr)[,"Seg_Filter"] <- "PASS"
	mcols (fit_gr)[,"mode"] <- NA_real_
	mcols (fit_gr)[fit_gr_data_gr.hits_data[,"queryHits"],"mode"] <- fit_gr_data_gr.hits_data[,"mode"]

	return (fit_gr)
}

intervals <- strsplit ("${intervals}", ",", fixed=T)[[1]]

interval_file <- gzfile ("${interval_bed}", 'rt')
data_interval <- read.table (file=interval_file,sep="\\t",header=F,stringsAsFactors=F)
names (data_interval) <- c("Chrom", "Chrom.Start", "Chrom.End")
head (data_interval)

data <- read.table (file="${loh_tsv}",sep="\\t",header=T,stringsAsFactors=F)
head (data)

seq_info <- GenomeInfoDb::Seqinfo (seqnames=data_interval %>% pull (Chrom),seqlengths=data_interval %>% pull (Chrom.End),genome="GRCh38.p12")
print (seq_info)

alt_haplotype_gr <- read.table (file="${alt_haplotype}",sep="\\t",header=T,stringsAsFactors=F) %>%
	makeGRangesFromDataFrame (seqinfo=seq_info)

centromere_gr <- read.table (file="${centromere}",sep="\\t",header=T,stringsAsFactors=F) %>%
	makeGRangesFromDataFrame (seqinfo=seq_info)

mappability_gr <- read.table (file="${mappability}",sep="\\t",header=T,stringsAsFactors=F) %>%
	makeGRangesFromDataFrame (seqinfo=seq_info,keep.extra.columns=T)

data_gr <- data %>%
	dplyr::mutate (LOH_Filter="PASS") %>%
	dplyr::select (Chrom,Pos,Plot_Freq,LOH_Filter) %>%
	makeGRangesFromDataFrame (seqinfo=seq_info,start.field="Pos",end.field="Pos",keep.extra.columns=T) %>%
	GenomeInfoDb::keepSeqlevels (intervals,pruning.mode="tidy")

print (data_gr)

# filter variants coming from "bad regions"
filter_gr <- IRanges::append (centromere_gr, alt_haplotype_gr)
data_gr_filter_gr.hits <- GenomicRanges::findOverlaps (data_gr, filter_gr)
mcols (data_gr)[queryHits (data_gr_filter_gr.hits),"LOH_Filter"] <- "bad_region"

data_gr_mappability_gr.hits <- GenomicRanges::findOverlaps (data_gr, mappability_gr)
data_gr_mappability_gr.hits_data <- as.data.frame (data_gr_mappability_gr.hits) %>%
	dplyr::mutate (score=!!mcols (mappability_gr)[subjectHits (data_gr_mappability_gr.hits),"score"]) %>%
	dplyr::group_by (queryHits) %>%
	dplyr::summarise (score=mean(score)) %>%
	data.frame

mcols (data_gr)[,"score"] <- NA_real_
mcols (data_gr)[data_gr_mappability_gr.hits_data[,"queryHits"],"score"] <- data_gr_mappability_gr.hits_data[,"score"]
print (data_gr)

mcols (data_gr)[which (mcols(data_gr)[,"LOH_Filter"] == "PASS" & (is.na (mcols (data_gr)[,"score"]) | mcols (data_gr)[,"score"] < 0.5)),"LOH_Filter"] <- "mappability"
print (data_gr)

seg_one_gr <- getLOHSegments (data_gr[which (mcols (data_gr)[,"LOH_Filter"] == "PASS")])
print ("seg_one_gr")
print (seg_one_gr)

# very small focal segments: remove all overlapping snps and rerun segmentation
minSegSize <- 10000000
data_gr_min_seg_gr.hits <- GenomicRanges::findOverlaps (data_gr,seg_one_gr[width (seg_one_gr)<minSegSize])
mcols (data_gr)[queryHits (data_gr_min_seg_gr.hits),"LOH_Filter"] <- "min_seg_size"

seg_two_gr <- getLOHSegments (data_gr[which (mcols (data_gr)[,"LOH_Filter"] == "PASS")])

print ("seg_two_gr")
print (seg_two_gr)

minSNPs <- 100

mcols (seg_two_gr)[which (mcols(seg_two_gr)[,"nbrOfLoci"]<minSNPs), "Seg_Filter"] <- "n_snps"

# merge similar neighboring segments

# now finally remove all segments with a very low number of SNPs per MB
minSNPsMB <- 100
mcols (seg_two_gr)[which (mcols (seg_two_gr)[,"Seg_Filter"] == "PASS" & ((mcols (seg_two_gr)[,"nbrOfLoci"] / width (seg_two_gr) * 1000000) < minSNPsMB)),"Seg_Filter"] <- "n_snps_mb"

write.table (as.data.frame (data_gr) %>% dplyr::rename (chromosome=seqnames,pos=start) %>% dplyr::select (chromosome,pos,Plot_Freq,LOH_Filter) %>% data.frame,file="${meta.sampleName}.LOH.snp_filter.tsv",sep="\\t",row.names=F,quote=F)
write.table (as.data.frame (seg_two_gr),file="${meta.sampleName}.LOH.segs.tsv",sep="\\t",row.names=F,quote=F)



	"""
}

process loh_seg_plot {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/LOH", mode: "copy"

	input:
		val (genome_build)
		val (intervals)
		tuple path (interval_bed), path (interval_bed_index)
		tuple val (meta), path (loh_snp), path (log_seg)

	script:
	"""#!/usr/bin/env Rscript

library (dplyr)
library (ggnewscale)
library (ggplot2)
library (gridExtra)

intervals <- strsplit ("${intervals}", ",", fixed=T)[[1]]

interval_file <- gzfile ("${interval_bed}", 'rt')
data_interval <- read.table (file=interval_file,sep="\t",header=F,stringsAsFactors=F)
names (data_interval) <- c("Chrom", "Start", "End")
head (data_interval)

data <- read.table (file="${log_snp}",sep="\t",header=T,stringsAsFactors=F)
head (data)

data_seg <- read.table (file="${loh_seg}",sep="\t",header=T,stringsAsFactors=F)
head (data_seg)

data_interval_plot <- data_interval %>%
	dplyr::mutate (Chrom=as.character (Chrom),End=as.numeric (End)) %>%
	dplyr::mutate (CumulativeStart=cumsum (End)-End) %>%
	dplyr::mutate (CumulativeMidpoint=CumulativeStart+((Start+End)/2)) %>%
	data.frame

data_plot <- data %>%
	dplyr::mutate (chromosome=as.character (chromosome)) %>%
	dplyr::filter (!is.na (Plot_Freq)) %>%
	dplyr::mutate (LOH_Filter=factor (LOH_Filter)) %>%
	dplyr::inner_join (data_interval_plot,by=c("chromosome"="Chrom"),suffix=c("",".Chrom")) %>%
	dplyr::mutate (Pos.Genome=pos+CumulativeStart) %>%
	data.frame

data_segments_plot <- data_seg %>%
	dplyr::mutate (chromosome=as.character (seqnames),Seg_Filter=factor (Seg_Filter)) %>%
	dplyr::mutate (upSeg=(0.5-mode)+0.5) %>%
	dplyr::mutate (lowSeg=1-upSeg) %>%
	dplyr::inner_join (data_interval_plot,by=c("chromosome"="Chrom"),suffix=c("",".Chrom")) %>%
	dplyr::mutate (Start.Genome=start+CumulativeStart,End.Genome=end+CumulativeStart) %>%
	data.frame

pdf (file="${meta.sampleName}.LOH.segments.genome.pdf",width=16,height=4)

ggplot (data_plot %>% filter (LOH_Filter=="PASS") %>% data.frame) +
	geom_point (aes (x=Pos.Genome,y=Plot_Freq),shape=".",colour="#ff80c3") +
	geom_segment (data=data_segments_plot %>% dplyr::filter (Seg_Filter=="PASS") %>% data.frame,aes(x=Start.Genome,y=upSeg,xend=End.Genome,yend=upSeg),colour="red") +
	geom_segment (data=data_segments_plot %>% dplyr::filter (Seg_Filter=="PASS") %>% data.frame,aes(x=Start.Genome,y=lowSeg,xend=End.Genome,yend=lowSeg),colour="red") +
	geom_hline (yintercept=0.5, color="darkgrey", alpha=0.5) +
	geom_vline (data=data_interval_plot,aes(xintercept=CumulativeStart)) +
	geom_text (data=data_interval_plot,aes(x=CumulativeMidpoint,y=1.1,label=Chrom),size=2) +
	coord_cartesian (ylim=c(0,1),expand=F,clip="off") +
	theme_bw () +
	theme (
		panel.grid.major.x=element_blank (),
		panel.grid.minor.x=element_blank (),
		axis.ticks.x=element_blank (),
		axis.text.x=element_blank (),
		plot.margin = unit(c(1,0.5,0.5,0.5), "cm")
	)

dev.off ()

plot_list <- vector ("list",length (intervals))


for ( i in seq_along (intervals) )
{
	plot_list[[i]] <- ggplot (data_plot %>% filter (chromosome==!!intervals[[i]],LOH_Filter=="PASS") %>% data.frame) +
		geom_point (aes(x=pos,y=Plot_Freq),shape=".",colour="#ff80c3") +
		geom_segment (data=data_segments_plot %>% filter (chromosome==!!intervals[[i]],Seg_Filter=="PASS") %>% data.frame,aes(x=start,y=upSeg,xend=end,yend=upSeg),colour="red") +
		geom_segment (data=data_segments_plot %>% filter (chromosome==!!intervals[[i]],Seg_Filter=="PASS") %>% data.frame,aes(x=start,y=lowSeg,xend=end,yend=lowSeg),colour="red") +
		ylim (0,1) +
		labs (title=chromosomes[[i]]) +
		theme_bw () +
		theme (
			panel.grid.major.x=element_blank (),
			panel.grid.minor.x=element_blank (),
			axis.text.x=element_blank ()

		)
}

# Need to do this outside of pdf call to prevent blank first page
p <- marrangeGrob (plot_list,nrow=1,ncol=1)

pdf (file="${meta.sampleName}.LOH.segments.chromosomes.pdf",width=9)
p
dev.off ()

# Plot again and include filtered snps and segments
loh_filter_colours <- setNames (rep ("#636363",length (levels (data_plot %>% pull (LOH_Filter)))), levels (data_plot %>% pull (LOH_Filter)))
loh_filter_colours[which (names (loh_filter_colours)=="PASS")] <- "#ff80c3"
print ("loh_filter_colours")
print (loh_filter_colours)
seg_filter_colours <- setNames (rep ("#636363",length (levels (data_segments_plot %>% pull (Seg_Filter)))),levels (data_segments_plot %>% pull (Seg_Filter)))
seg_filter_colours[which (names (seg_filter_colours)=="PASS")] <- "red"

pdf (file="${meta.sampleName}.LOH.segments.genome.full.pdf",width=16,height=4)

ggplot (data_plot) +
	geom_point (aes (x=Pos.Genome,y=Plot_Freq,colour=LOH_Filter),shape=".",show.legend=F) +
	scale_colour_manual (values=loh_filter_colours) +
	ggnewscale::new_scale_color () +
	geom_segment (data=data_segments_plot,aes(x=Start.Genome,y=upSeg,xend=End.Genome,yend=upSeg,colour=Seg_Filter),show.legend=F) +
	geom_segment (data=data_segments_plot,aes(x=Start.Genome,y=lowSeg,xend=End.Genome,yend=lowSeg,colour=Seg_Filter),show.legend=F) +
	scale_colour_manual (values=seg_filter_colours) +
	geom_hline (yintercept=0.5, color="darkgrey", alpha=0.5) +
	geom_vline (data=data_interval_plot,aes(xintercept=CumulativeStart)) +
	geom_text (data=data_interval_plot,aes(x=CumulativeMidpoint,y=1.1,label=Chrom),size=2) +	
	coord_cartesian (ylim=c(0,1),expand=F,clip="off") +
	theme_bw () +
	theme (
		panel.grid.major.x=element_blank (),
		panel.grid.minor.x=element_blank (),
		axis.ticks.x=element_blank (),
		axis.text.x=element_blank (),
		plot.margin = unit(c(1,0.5,0.5,0.5), "cm")
	)

dev.off ()

plot_list <- vector ("list",length (intervals))


for ( i in seq_along (intervals) )
{
	plot_list[[i]] <- ggplot (data_plot %>% filter (chromosome==!!intervals[[i]]) %>% data.frame) +
		geom_point (aes(x=pos,y=Plot_Freq,colour=LOH_Filter),shape=".",show.legend=F) +
		scale_colour_manual (values=loh_filter_colours) +
		ggnewscale::new_scale_color () +
		geom_segment (data=data_segments_plot %>% filter (chromosome==!!intervals[[i]]) %>% data.frame,aes(x=start,y=upSeg,xend=end,yend=upSeg,colour=Seg_Filter)) +
		geom_segment (data=data_segments_plot %>% filter (chromosome==!!intervals[[i]]) %>% data.frame,aes(x=start,y=lowSeg,xend=end,yend=lowSeg,colour=Seg_Filter)) +
		ylim (0,1) +
		labs (title=intervals[[i]]) +
		theme_bw () +
		theme (
			panel.grid.major.x=element_blank (),
			panel.grid.minor.x=element_blank (),
			axis.text.x=element_blank ()

		)
}

# Need to do this outside of pdf call to prevent blank first page
p <- marrangeGrob (plot_list,nrow=1,ncol=1)

pdf (file="${meta.sampleName}.LOH.segments.chromosomes.full.pdf",width=9)
p
dev.off ()


	"""
}


