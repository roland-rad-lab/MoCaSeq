
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
	dplyr::mutate (Chrom=as.character (Chrom), End=as.numeric (End)) %>%
	dplyr::mutate (CumulativeStart=cumsum (End)-End) %>%
	dplyr::mutate (CumulativeEnd=cumsum (End)) %>%
	dplyr::mutate (CumulativeMidpoint=(CumulativeStart+CumulativeEnd)/2) %>%
	data.frame

data_plot <- data %>%
	dplyr::inner_join (data_interval_plot,by="Chrom",suffix=c("",".Chrom")) %>%
	dplyr::mutate (Pos.Genome=Pos+CumulativeStart) %>%
	data.frame

plot_types <- setNames (c("Plot_Freq","Tumor_Freq","Normal_Freq"),c("adjusted","raw","germline"))

for ( i in seq_along (plot_types) )
{
	g <- ggplot (data_plot) +
		geom_point (aes_string(x="Pos.Genome",y=plot_types[i]),shape=".",colour="#ffa8d7") +
		geom_vline (data=data_interval_plot,aes(xintercept=CumulativeStart)) +
		geom_text (data=data_interval_plot,aes(x=CumulativeMidpoint,y=1.05,label=Chrom),size=2) +
		coord_cartesian (xlim=c(0,data_interval_plot %>% pull (CumulativeEnd) %>% max ()),ylim=c(0,1),expand=F,clip="off") +
		xlab ("Genome") +
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
			geom_point (aes_string(x="Pos",y=plot_types[i]),shape=".",colour="#ffa8d7") +
			scale_x_continuous (labels=scales::number_format (big.mark=",",scale=1e-06,suffix=" Mb",accuracy=0.1)) +
			coord_cartesian (xlim=c(0,data_interval_plot %>% filter (Chrom==!!chromosomes[[j]]) %>% pull (End)),ylim=c(0,1)) +
			labs (title=chromosomes[[j]]) +
			theme_bw () +
			theme (
				panel.grid.major.x=element_blank (),
				panel.grid.minor.x=element_blank (),
				plot.margin=unit (c(5.5,25.5,5.5,5.5),"pt")
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
touch ${meta.sampleName}.LOH.germline.chromosomes.pdf
touch ${meta.sampleName}.LOH.germline.genome.pdf
touch ${meta.sampleName}.LOH.raw.chromosomes.pdf
touch ${meta.sampleName}.LOH.raw.genome.pdf

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
		tuple val (meta), path ("${meta.sampleName}_LOH_SNPs.tsv"), path ("${meta.sampleName}_LOH_Segments.tsv"), emit: result

	script:
	"""#!/usr/bin/env Rscript
library (data.table)
library (dplyr)
library (GenomicRanges)
library (PSCBS)
library (zoo)

# Following method by NikdAK

ConvertGenomicCords = function(dat,chrom.sizes,Start,CopyNumber,Chromosome)
{
  Outlist = list()
  cn <- data.frame()
  borders <- c()
  ChromBorders = c()
  FirstPosition = c()
  last = 0
  for(i in names(chrom.sizes))
  {
    cur <- dat[dat[,Chromosome]==i,c(Start,CopyNumber)]
    cn <- rbind(cn,data.frame(Chromosome = i,position=cur[,Start]+last,copy=cur[,CopyNumber]))
    borders <- c(borders,last)
    last = last + chrom.sizes[i]
    ChromBorders = c(ChromBorders,last)
    FirstPosition = c(FirstPosition,min(dat[dat[,Chromosome]==i,Start]))
  }
  names(ChromBorders) = names(chrom.sizes)
  Outlist[["CN"]]=cn
  Outlist[["ChromosomeBorders"]] = ChromBorders
  Outlist[["FirstChromPosition"]] = FirstPosition
  Outlist[["NonProcessed"]] = dat
  return(Outlist)
}

ProcessCountData2 = function(countdata="",chrom.sizes=chrom.sizes,method=""){
  Outlist = list()
  countdata = fread(countdata,header=T,sep="\t")
  SetVariableNames(method)
  countdata <- data.frame(countdata)
  Outlist = ConvertGenomicCords(countdata,chrom.sizes,Start,CopyNumber,Chromosome)
  return(Outlist)
}

MirrorAndDuplicate <- function(lohDT){
  # first we will "mirror&duplicate" all points:
  # - points can be categorized in up or below 0.5
  # - now all values below 0.5 will be copied to the upper but keeping the same distance to 0.5 as it was before (mirroring)
  # - the same happens for the above group so in the end the data points will be doubled

  upperDT <- copy(lohDT)
  upperDT[, LOHscore := Plot_Freq]
  upperDT[Plot_Freq <= 0.5, LOHscore := 0.5+(0.5-Plot_Freq)] # move all values from below 0.5 up
  lowerDT <- copy(upperDT) # now copy the entire set back down (duplicating the data points)
  lowerDT[, LOHscore := 1-LOHscore]

  dupDT <- rbind(upperDT,lowerDT) # the duplicated lohDT
  return(dupDT)
}

AssignSNPids <- function(lohSegs,lohDT){
  lohSegsGR <- makeGRangesFromDataFrame(lohSegs[, .(chr=Chrom, start,end)])
  lohGR <- makeGRangesFromDataFrame(lohDT[, .(chr=Chrom, start=Pos, end=Pos)])
  hits <- findOverlaps(lohSegsGR, lohGR)
  lohDT[subjectHits(hits), segID := lohSegs[queryHits(hits), segID]]
  lohDT[, snpID := .I] # also SNP ID

  return(lohDT)
}


# CBS segmentation algorithm
GetSegments <- function(dupDT){
  # run segmentation
  data <- data.frame(dupDT[LOHscore >= 0.5, .(chromosome=as.character(Chrom), x=Pos, y=LOHscore)])
  data <- data[, c("chromosome", "x", "y")]
  rownames(data) <- paste0("SEG", 1:nrow(data))
  gaps <- findLargeGaps(data, minLength = 2e+06) # smallest centromere is 2100000
  knownSegments <- gapsToSegments(gaps)
  fit <- segmentByCBS(data, knownSegments = knownSegments, verbose = -10, seed=T)

  lohSegs <- as.data.table(fit)
  setnames(lohSegs, "chromosome", "Chrom")
  lohSegs[, Chrom := as.character(Chrom)]
  lohSegs <- lohSegs[!is.na(mean)]
  lohSegs[, width := round(end-start)]
  lohSegs[, segID := sequence(.N), by = Chrom] # segment ID
  lohSegs[, mid := start+((end-start)/2), by=segID] # center for plot labels

  # assign (mirror&duplicated) SNPs to segments
  lohSegsGR <- makeGRangesFromDataFrame(lohSegs[, .(chr=Chrom, start,end)])
  dupGR <- makeGRangesFromDataFrame(dupDT[, .(chr=Chrom, start=Pos, end=Pos)])
  hits <- findOverlaps(lohSegsGR, dupGR)
  dupDT[subjectHits(hits), segID := lohSegs[queryHits(hits), segID]]

  # calculate the mode value for this segment
  modeVals <- dupDT[, density(LOHscore)\$x[which.max(density(LOHscore)\$y)], by=.(Chrom, segID)]
  names(modeVals) <- c("Chrom","segID", "mode")
  #modeVals[, Chrom := as.integer(Chrom)]

  lohSegs <- merge(lohSegs, modeVals, by=c("Chrom", "segID")) # assign the calculcated mode values for each segment

  lohSegs[, mean := NULL] # this is meaningless

  # assign the new mode segment values to the upper and lower region
  lohSegs[, upSeg := (0.5-mode)+0.5]
  lohSegs[, lowSeg := 1-upSeg]
  lohSegs[lowSeg > upSeg, c("upSeg", "lowSeg") := .(lowSeg, upSeg)]

  return(lohSegs)
}

#Get_LOH_segments <- function(name,seqType,species,filterBAF=T,denovoBAF=F,minLOHsegdif=0.03,customData=NULL,MinDifToMerge=0.05,minSegSize=100000)
Get_LOH_segments <- function(name,seqType,species,data,data_gr,chrom.sizes,filterBAF=T,denovoBAF=F,minLOHsegdif=0.03,customData=NULL,MinDifToMerge=0.05,minSegSize=100000)
{

  # denovoBAF: take SNP BAFs from table or recalculate
  # filterBAF: apply snp filtering (centromeres, haplotypes, mappability) --> currently only for human
  # minLOHsegdif: up/low segment with differences smaller than this will be set to 0.5
  # MinDifToMerge: 0.05 # distance (BAF difference) used to merge neighbouring segments
  # minSegSize <- 100000 # minimal length of called segments

  # minimal number of total SNPs and SNPs per megabase
  if(seqType == "WES"){
    minSNPsMB <- 3
    minSNPs <- 5
  } else if(seqType == "WGS") {
    minSNPs <- 100
    minSNPsMB <- 100
  }

####
#  chrom.sizes = DefineChromSizes(species)
#  if (species=="human"){
#    chromosomes=21
#    nXbreaks <- 10 # number of breaks on x-axis
#
#    # data calculated in: "/home/rad/Documents/small_scripts/MoCaSeq_DevCode/Cohort_GenerateOverlay_FilterNoise.R"
#    load("/media/rad/SSD2/MoCaSeq_runs/AGRad_test/loh_seg/hg38_cents-haplo-mappability.RData")
#  } else if (species=="Mouse"){
#    chromosomes=19
#    nXbreaks <- 10 # number of breaks on x-axis
#    load("/media/rad/HDD1/Lookups/mm10_cents-haplo-mappability.RData")
#  }
#
#
#  if(!denovoBAF){
#    # OLD
#    data=paste0(name,"/results/LOH/",name,".VariantsForLOH.txt")
#    LOHDat = ProcessCountData2(data,chrom.sizes,"LOH")
#    lohDT <- data.table(LOHDat[[4]])
#    lohDT <- lohDT[!Chrom %in% c("X", "Y")]
#    lohDT[, Chrom := as.numeric(Chrom)]
#  } else {
#
#    # NEW, this takes a long time
#    tumorBAM <- paste(name,"/results/bam/",name,".Tumor.bam", sep="")
#    tumor <- fread(paste(name,"/results/Mutect2/",name,".Tumor.Mutect2.Positions.txt", sep=""), header=T, sep="\t",
#                   col.names=c("Chrom", "Pos", "Ref", "Alt", "Tumor_Freq", "Tumor_RefCount", "Tumor_AltCount", "MapQ", "BaseQ"))
#    normal <- fread(paste(name,"/results/Mutect2/",name,".Normal.Mutect2.Positions.txt", sep=""), header=T, sep="\t",
#                    col.names=c("Chrom", "Pos", "Ref", "Alt", "Normal_Freq", "Normal_RefCount", "Normal_AltCount", "MapQ", "BaseQ"))
#
#    # filter normal (keep tumor as it is)
#    normal <- normal[Chrom %in% c(1:22)]
#    normal <- normal[nchar(paste0(Ref, Alt)) == 2]
#    normal <- normal[MapQ >= 60]
#    normal <- normal[(Normal_RefCount + Normal_AltCount) >= 10]
#    normal <- normal[Normal_Freq %between% c(0.3, 0.7)]
#    normal[, id := paste0(Chrom, ":", Pos, "-", Pos)]
#    normal[, UniquePos := paste(Chrom, Pos, Pos, sep="_")]
#    tumor[, UniquePos := paste(Chrom, Pos, Pos, sep="_")]
#
#    # Scan the tumor bam file for coverage
#    library(Rsamtools)
#    which <- GRanges(seqnames = normal\$Chrom, ranges = IRanges(normal\$Pos,normal\$Pos))
#    what <- c("mapq")
#    param <- ScanBamParam(which = which, what = what)
#    tumorList <- scanBam(tumorBAM, param=param)
#    tumorList <- lapply(tumorList, unlist)
#    tumorList <- lapply(tumorList, function(x) length(x[x >= 60])) # apply MAPQ filtering
#    tumorRaw <- data.table(TumorCounts=tumorList)
#    tumorRaw[, id := names(tumorList)]
#
#    variants <- merge(normal, tumorRaw, by="id", all.x=T)
#    variants <- merge(variants, tumor[, .(UniquePos, Tumor_Freq)], all.x=T, by="UniquePos")
#    variants[, Chrom := as.numeric(Chrom)]
#    #saveRDS(variants, "/home/rad/Downloads/temp6/DS36_variants_unfiltered.rds")
#    #variants <- readRDS("/home/rad/Downloads/temp6/DS36_variants_unfiltered.rds")
#    variants <- variants[TumorCounts >= 10]
#    variants[is.na(Tumor_Freq), Tumor_Freq := 0]
#
#    #ggplot(variants[Chrom == 6], aes(Pos, Tumor_Freq)) + geom_point()
#    variants[, Dif_Freq := Tumor_Freq-Normal_Freq]
#    variants[, Tumor_Freq := Dif_Freq]
#    lohDT <- copy(variants)
#  }
#
####

  lohDT = data.table (data)

  # reset the "mouse correction procedure" for human
  if ( species=="human" )
  {
    lohDT[, Plot_Freq := Tumor_Freq]
  }

  # generate a copy with all SNPs (and merge it back to keep removed SNPs but with a lighter color)
  lohDT.all <- copy(lohDT)

  lohGR <- makeGRangesFromDataFrame(lohDT[, .(chr=Chrom, start=Pos, end=Pos)])

  # filter SNPS for centromeres, low mappability regions and alternative haplotypes
  if(filterBAF){
    # filter variants coming from "bad regions"

####
#    filter.gr <- append(cents.gr, haplo.gr)
#    hits <- findOverlaps(lohGR, filter.gr)
#    lohDT[queryHits(hits), remove := "yes"]
#
#    hits <- findOverlaps(lohGR, mappab.gr)
#    BAFmappa <- cbind(lohDT[queryHits(hits)], mappab[subjectHits(hits), score])
#    BAFmappa <- BAFmappa[, mean(V2), by=UniquePos]
#    setnames(BAFmappa, "V1", "mappability")
#    lohDT <- merge(lohDT, BAFmappa, by="UniquePos", all.x=T, sort=F)
#    lohDT[is.na(mappability), mappability := 0]
#
#    lohDT <- lohDT[is.na(remove)]
#    lohDT <- lohDT[mappability > 0.5]
####
    lohDT <- lohDT[mcols (data_gr)[mcols(data_gr)[,"LOH_Filter"]=="PASS" & !is.na (mcols(data_gr)[,"score"]) & mcols(data_gr)[,"score"] > 0.5,"row_index"]]
    lohGR <- makeGRangesFromDataFrame(lohDT[, .(chr=Chrom, start=Pos, end=Pos)])
  }

  dupDT <- MirrorAndDuplicate(lohDT)
  lohSegs <- GetSegments(dupDT)

  # if segments and dots are removed at the start/end of the chromosome, the points in the middle would shift and look like the start of the chromosome
  # here we generate dummy data to keep the lengths as they were before
  startDummy <- data.table(Chrom = names(chrom.sizes), Pos = 0, Plot_Freq=0)
  endDummy <- data.table(Chrom = names(chrom.sizes), Pos = chrom.sizes, Plot_Freq=0)
  chromDummy <- rbind(startDummy, endDummy)

  # very small focal segments: remove all overlapping SNPs and rerun segmentation
  removeSegs <- lohSegs[width < minSegSize]
  if(nrow(removeSegs) > 0){
    removeSegs.gr <- makeGRangesFromDataFrame(removeSegs[, .(chr=Chrom, start, end)])
    hits <- findOverlaps(lohGR, removeSegs.gr)
    lohDT <- lohDT[!queryHits(hits)]

    dupDT <- MirrorAndDuplicate(lohDT)
  }

  lohSegs <- GetSegments(dupDT)

  if(nrow(lohSegs) == 0){print("No LOH segments found at checkpoint A (change the parameter \\"minSegSize\\")");quit ('no')}

  # re-assign segID to the SNPs
  lohDT <- AssignSNPids(lohSegs, lohDT)

  # remove segments supported only by few SNP, as well as the corresponding SNPs
  snpcountDT <- lohDT[, .N, by=.(Chrom, segID)]
  lohSegs <- merge(lohSegs, snpcountDT[, .(Chrom, segID, nSNPs=N)], by=c("Chrom", "segID")) # also add it to segments
  removeSegs <- snpcountDT[N <= minSNPs, paste(Chrom, segID, sep="+")]
  removeSNPs <- lohDT[paste(Chrom, segID,sep="+") %in% removeSegs, snpID]

  lohSegs <- lohSegs[!paste(Chrom, segID, sep="+") %in% removeSegs]
  lohDT <- lohDT[!snpID %in% removeSNPs]

  if(nrow(lohSegs) == 0){print("No LOH segments found at checkpoint B (change the parameter \\"minSNPs\\")");quit ('no')}

  # merge similar neighboring segments
  lohSegs[, DifNextSeg := abs(upSeg - data.table::shift(upSeg)), by = Chrom]
  lohSegs[, DistNextSeg := abs(start - data.table::shift(end)), by = Chrom]
  lohSegs[, mergeID := as.numeric(.I)]
  lohSegs[DifNextSeg < MinDifToMerge & DistNextSeg < 10000000, mergeID := NA]
  lohSegs[, mergeID := na.locf(mergeID)]

  lohSegs <- unique(lohSegs[, .(Chrom, segID=as.integer(mean(segID)),
                                start=min(start), end=max(end), mid=mean(mid), mode=mean(mode),
                                upSeg=mean(upSeg), lowSeg=mean(lowSeg), width=sum(width), nSNPs=sum(nSNPs)), by=mergeID])
  lohSegs[, mergeID := NULL]

  # reassign segment IDs
  lohSegs[, segID := sequence(.N), by = Chrom]

  # re-assign segID to the SNPs
  lohDT <- AssignSNPids(lohSegs, lohDT)

  # now finally remove all segments with a very low number of SNPs per MB
  lohSegs[, nSNPperMB := nSNPs / width * 1000000]
  removeSegs <- lohSegs[nSNPperMB < minSNPsMB, paste(Chrom, segID, sep="+")]
  removeSNPs <- lohDT[paste(Chrom, segID,sep="+") %in% removeSegs, snpID]
  lohSegs <- lohSegs[!paste(Chrom, segID, sep="+") %in% removeSegs]
  lohDT <- lohDT[!snpID %in% removeSNPs]

  if(nrow(lohSegs) == 0){print("No LOH segments found at checkpoint C (change the parameter \\"minSNPsMB\\")");quit ('no')}

  # reassign segment IDs
  lohSegs[, segID := sequence(.N), by = Chrom]

  # re-assign segID to the SNPs
  lohDT <- AssignSNPids(lohSegs, lohDT)

  # add removed SNPs back but mark them
  lohDT.removed <- lohDT.all[!UniquePos %in%lohDT\$UniquePos]

  lohSegs[, segDif := upSeg-lowSeg]
  lohSegs[segDif < minLOHsegdif, `:=` (upSeg=0.5, lowSeg=0.5)]

  return(list(lohDT=lohDT, lohSegs=lohSegs))
}

intervals <- strsplit ("${intervals}", ",", fixed=T)[[1]]

interval_file <- gzfile ("${interval_bed}", 'rt')
data_interval <- read.table (file=interval_file,sep="\\t",header=F,stringsAsFactors=F) %>%
	dplyr::rename (Chrom=V1,Chrom.Start=V2,Chrom.End=V3) %>%
	dplyr::mutate (Chrom=as.character (Chrom)) %>%
	data.frame
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
	dplyr::mutate (row_index=row_number (),LOH_Filter="PASS") %>%
	dplyr::select (Chrom,Pos,Tumor_Freq,Plot_Freq,LOH_Filter,row_index) %>%
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

chrom.sizes <- setNames (data_interval %>% dplyr::filter (Chrom %in% intervals) %>% pull (Chrom.End), data_interval %>% filter (Chrom %in% intervals) %>% pull (Chrom))

result <- Get_LOH_segments ("${meta.sampleName}",toupper ("${meta.seqType}"), "${meta.organism}",data,data_gr,chrom.sizes)

write.table (result[["lohDT"]],file="${meta.sampleName}_LOH_SNPs.tsv",sep="\\t",row.names=F,quote=F)
write.table (result[["lohSegs"]],file="${meta.sampleName}_LOH_Segments.tsv",sep="\\t",row.names=F,quote=F)

	"""

	stub:
	"""#!/usr/bin/env bash

if [[ "${params.stub_json_map?.loh_matched_segment}" == "null" ]]; then
        cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/LOH/${meta.sampleName}_LOH_SNPs.tsv .
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/LOH/${meta.sampleName}_LOH_Segments.tsv .
fi

touch ${meta.sampleName}_LOH_SNPs.tsv
touch ${meta.sampleName}_LOH_Segments.tsv

	"""
}

process loh_seg_plot {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/LOH", mode: "copy"

	input:
		val (genome_build)
		val (intervals)
		tuple path (interval_bed), path (interval_bed_index)
		tuple val (meta), path (loh_snp), path (loh_seg)

	output:
		path ("${meta.sampleName}.*.pdf")

	script:
	"""#!/usr/bin/env Rscript

library (dplyr)
library (ggnewscale)
library (ggplot2)
library (ggrepel)
library (gridExtra)

intervals <- strsplit ("${intervals}", ",", fixed=T)[[1]]

interval_file <- gzfile ("${interval_bed}", 'rt')
data_interval <- read.table (file=interval_file,sep="\t",header=F,stringsAsFactors=F)
names (data_interval) <- c("Chrom", "Start", "End")
head (data_interval)

data <- read.table (file="${loh_snp}",sep="\t",header=T,stringsAsFactors=F)
head (data)

data_seg <- read.table (file="${loh_seg}",sep="\t",header=T,stringsAsFactors=F)
head (data_seg)

data_interval_plot <- data_interval %>%
	dplyr::mutate (Chrom=as.character (Chrom),End=as.numeric (End)) %>%
	dplyr::mutate (CumulativeStart=cumsum (End)-End) %>%
	dplyr::mutate (CumulativeEnd=cumsum (End)) %>%
	dplyr::mutate (CumulativeMidpoint=(CumulativeStart+CumulativeEnd)/2) %>%
	dplyr::filter (Chrom %in% !!intervals) %>%
	data.frame

data_plot <- data %>%
	dplyr::mutate (Chrom=as.character (Chrom)) %>%
	dplyr::filter (!is.na (Plot_Freq)) %>%
	dplyr::mutate (LOH_Filter="PASS") %>%
	dplyr::mutate (LOH_Filter=factor (LOH_Filter)) %>%
	dplyr::inner_join (data_interval_plot,by="Chrom",suffix=c("",".Chrom")) %>%
	dplyr::mutate (Pos.Genome=Pos+CumulativeStart) %>%
	data.frame

data_segments_plot <- data_seg %>%
	dplyr::mutate (Seg_Filter="PASS") %>%
	dplyr::mutate (Chrom=as.character (Chrom),Seg_Filter=factor (Seg_Filter)) %>%
	dplyr::inner_join (data_interval_plot,by="Chrom",suffix=c("",".Chrom")) %>%
	dplyr::mutate (Start.Genome=start+CumulativeStart,End.Genome=end+CumulativeStart) %>%
	data.frame

pdf (file="${meta.sampleName}.LOH.segments.genome.pdf",width=16,height=4)

ggplot (data_plot %>% filter (LOH_Filter=="PASS") %>% data.frame) +
	geom_point (aes (x=Pos.Genome,y=Plot_Freq),shape=".",colour="#ffa8d7") +
	geom_segment (data=data_segments_plot %>% dplyr::filter (Seg_Filter=="PASS") %>% data.frame,aes(x=Start.Genome,y=upSeg,xend=End.Genome,yend=upSeg),colour="red") +
	geom_segment (data=data_segments_plot %>% dplyr::filter (Seg_Filter=="PASS") %>% data.frame,aes(x=Start.Genome,y=lowSeg,xend=End.Genome,yend=lowSeg),colour="red") +
	geom_hline (yintercept=0.5, color="darkgrey", alpha=0.5) +
	geom_vline (data=data_interval_plot,aes(xintercept=CumulativeStart)) +
	geom_text (data=data_interval_plot,aes(x=CumulativeMidpoint,y=1.05,label=Chrom),size=2) +
	coord_cartesian (xlim=c(0,data_interval_plot %>% pull (CumulativeEnd) %>% max ()),ylim=c(0,1),expand=F,clip="off") +
	xlab ("Genome") +
	ylab ("BAF") +
	theme_bw () +
	theme (
		panel.grid.major.x=element_blank (),
		panel.grid.minor.x=element_blank (),
		axis.ticks.x=element_blank (),
		axis.text.x=element_blank (),
		plot.margin=unit (c(1,0.5,0.5,0.5), "cm")
	)

dev.off ()

plot_list <- vector ("list",length (intervals))


for ( i in seq_along (intervals) )
{
	plot_list[[i]] <- ggplot (data_plot %>% filter (Chrom==!!intervals[[i]],LOH_Filter=="PASS") %>% data.frame) +
		geom_point (aes(x=Pos,y=Plot_Freq),shape=".",colour="#ffa8d7") +
		geom_segment (data=data_segments_plot %>% filter (Chrom==!!intervals[[i]],Seg_Filter=="PASS") %>% data.frame,aes(x=start,y=upSeg,xend=end,yend=upSeg),colour="red") +
		geom_segment (data=data_segments_plot %>% filter (Chrom==!!intervals[[i]],Seg_Filter=="PASS") %>% data.frame,aes(x=start,y=lowSeg,xend=end,yend=lowSeg),colour="red") +
		scale_x_continuous (labels=scales::number_format (big.mark=",",scale=1e-06,suffix=" Mb",accuracy=0.1)) +
		coord_cartesian (xlim=c(0,data_interval_plot %>% filter (Chrom==!!intervals[[i]]) %>% pull (End)),ylim=c(0,1)) +
		labs (title=intervals[[i]]) +
		ylab ("BAF") +
		theme_bw () +
		theme (
			panel.grid.major.x=element_blank (),
			panel.grid.minor.x=element_blank (),
			plot.margin=unit (c(5.5,25.5,5.5,5.5),"pt")
		)
}

# Need to do this outside of pdf call to prevent blank first page
p <- marrangeGrob (plot_list,nrow=1,ncol=1)

pdf (file="${meta.sampleName}.LOH.segments.chromosomes.pdf",width=9)
p
dev.off ()

# This time we include labels for interesting segments
plot_list <- vector ("list",length (intervals))


for ( i in seq_along (intervals) )
{
	plot_list[[i]] <- ggplot (data_plot %>% filter (Chrom==!!intervals[[i]],LOH_Filter=="PASS") %>% data.frame) +
		geom_point (aes(x=Pos,y=Plot_Freq),shape=".",colour="black") +
		geom_segment (data=data_segments_plot %>% filter (Chrom==!!intervals[[i]],Seg_Filter=="PASS") %>% data.frame,aes(x=start,y=upSeg,xend=end,yend=upSeg),colour="red") +
		geom_segment (data=data_segments_plot %>% filter (Chrom==!!intervals[[i]],Seg_Filter=="PASS") %>% data.frame,aes(x=start,y=lowSeg,xend=end,yend=lowSeg),colour="red") +
		geom_label_repel (data=data_segments_plot %>% filter (Chrom==!!intervals[[i]],Seg_Filter=="PASS",upSeg<0.8) %>% data.frame,aes(mid, upSeg, label=segID),nudge_y=0.1,size=3) +
		geom_label_repel (data=data_segments_plot %>% filter (Chrom==!!intervals[[i]],Seg_Filter=="PASS",upSeg>0.8) %>% data.frame,aes(mid, upSeg, label=segID),nudge_y=-0.1,size=3) +
		scale_x_continuous (labels=scales::number_format (big.mark=",",scale=1e-06,suffix=" Mb",accuracy=0.1)) +
		coord_cartesian (xlim=c(0,data_interval_plot %>% filter (Chrom==!!intervals[[i]]) %>% pull (End)),ylim=c(0,1)) +
		labs (title=intervals[[i]]) +
		ylab ("BAF") +
		theme_bw () +
		theme (
			panel.grid.major.x=element_blank (),
			panel.grid.minor.x=element_blank (),
			plot.margin=unit (c(5.5,25.5,5.5,5.5),"pt")
		)
}

# Need to do this outside of pdf call to prevent blank first page
p <- marrangeGrob (plot_list,nrow=1,ncol=1)

pdf (file="${meta.sampleName}.LOH.segments.chromosomes.labelled.pdf",width=9)
p
dev.off ()


# Plot again and include filtered snps and segments
loh_filter_colours <- setNames (rep ("#dedede",length (levels (data_plot %>% pull (LOH_Filter)))), levels (data_plot %>% pull (LOH_Filter)))
loh_filter_colours[which (names (loh_filter_colours)=="PASS")] <- "#ffedf7"

seg_filter_colours <- setNames (rep ("blue",length (levels (data_segments_plot %>% pull (Seg_Filter)))),levels (data_segments_plot %>% pull (Seg_Filter)))
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
	geom_text (data=data_interval_plot,aes(x=CumulativeMidpoint,y=1.05,label=Chrom),size=2) +
	coord_cartesian (xlim=c(0,data_interval_plot %>% pull (CumulativeEnd) %>% max ()),ylim=c(0,1),expand=F,clip="off") +
	theme_bw () +
	xlab ("Genome") +
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
	plot_list[[i]] <- ggplot (data_plot %>% filter (Chrom==!!intervals[[i]]) %>% data.frame) +
		geom_point (aes(x=Pos,y=Plot_Freq,colour=LOH_Filter),shape=".",show.legend=F) +
		scale_colour_manual (values=loh_filter_colours) +
		ggnewscale::new_scale_color () +
		geom_segment (data=data_segments_plot %>% filter (Chrom==!!intervals[[i]]) %>% data.frame,aes(x=start,y=upSeg,xend=end,yend=upSeg,colour=Seg_Filter),show.legend=F) +
		geom_segment (data=data_segments_plot %>% filter (Chrom==!!intervals[[i]]) %>% data.frame,aes(x=start,y=lowSeg,xend=end,yend=lowSeg,colour=Seg_Filter),show.legend=F) +
		scale_colour_manual (values=seg_filter_colours) +
		scale_x_continuous (labels=scales::number_format (big.mark=",",scale=1e-06,suffix=" Mb",accuracy=0.1)) +
		coord_cartesian (xlim=c(0,data_interval_plot %>% filter (Chrom==!!intervals[[i]]) %>% pull (End)),ylim=c(0,1)) +
		labs (title=intervals[[i]]) +
		theme_bw () +
		theme (
			panel.grid.major.x=element_blank (),
			panel.grid.minor.x=element_blank (),
			plot.margin=unit (c(5.5,25.5,5.5,5.5),"pt")
		)
}

# Need to do this outside of pdf call to prevent blank first page
p <- marrangeGrob (plot_list,nrow=1,ncol=1)

pdf (file="${meta.sampleName}.LOH.segments.chromosomes.full.pdf",width=9)
p
dev.off ()

	"""

	stub:
	"""#!/usr/bin/env bash

touch ${meta.sampleName}.LOH.segments.chromosomes.pdf
touch ${meta.sampleName}.LOH.segments.chromosomes.full.pdf
touch ${meta.sampleName}.LOH.segments.chromosomes.labelled.pdf
touch ${meta.sampleName}.LOH.segments.genome.pdf
touch ${meta.sampleName}.LOH.segments.genome.full.pdf
	"""
}


