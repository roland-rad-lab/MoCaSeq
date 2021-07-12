
process loh_matched {
	tag "${meta.sampleName}"

	input:
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

	publishDir "${params.output_base}/${meta.sampleName}/results/LOH", mode: "copy"

	input:
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

cp ${params.stub_dir}/${meta.sampleName}/results/LOH/${meta.sampleName}.VariantsForLOH.txt .

	"""

}



