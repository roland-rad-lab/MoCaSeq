#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// print help info TODO
if (params.help != null){
	println """\
MoCaSeq cancer genome sequencing pipeline for mouse and human
======================================================================================
Extensive analysis of cancer sequencing data from human or mouse featuring mapping,
 remapping, mutation calling, LOH analysis, 
Employs: cnv-kit, bwa_mem, mutect2, 

Usage: 
nextflow run roland-rad-lab/MoCaSeq -r human-pipeline-nextflow-2 --input [PARAMETERS] 

1. General Nextflow parameters
-entry				set specific workflow as entry point by name
-profile			use nextflow profiles defined in config files
-work-dir			set working directory

2. MoCaSeq specific parameters
2.1. I/O
--output_base			specify output path prefix
--input				set input .tsv file with sample information 
				 (see git repo for details)
--qc_dir			|

2.2. Config
--cache_base			directory that hold pipeline cache dirs
--script_base			directory holding some pipeline scripts inside containers.
				Make sure to bind MoCaSeq/repository/:/var/pipeline/repository
				or custom path with respective change in your mocaseq.config. 
--custom_config_base		directory to look for mocaseq.config, genomes.config
				and genome_annotations.config. These *.config files 
				need to be in a subdir '<custom_config_version>/pipeline'.
				For templates see /conf in git repo. May be https url.
--custom_config_version		config version found in <custom_config_base>
--test_config_genome_base	directory for test *.config files like --custom_config_base
				IMPORTANT: this overwrites --custom_config_base
--test_config_genome_version	test config version like --custom_config_version
				IMPORTANT: this overwrites --custom_config_version
--gatk.jar			Path to Genome Analysis Toolkit .jar file.
				Please use /opt/gatk-<version>/gatk.jar with version one of
				(4.2.0.0, 4.1.7.0, 4.1.3.0, 4.1.0.0) for this. Alternatively
				adapt the container for another version of gatk.
--tiny				skip pipeline steps that would fail with few data
		
2.3. Genomes
--genomes_base			|
--genomes			|
--genome_annotations		|
--genome_build.human		|
--genome_build.mouse		|
--stub_json			Optionally load json map to control the behaviour of 
				stubs (cp vs touch)

2.4. Panel of normals analysis
--pon_dir			|
--pon_sample			|
--pon_tsv			|
--pon_bed			|
--pon_dry			|
--pon_resolution		|

2.5. Tool parameters
--cnv_kit			ch_centre assignment
--cnv_kit.centre		|
--cnv_kit.centre.tokenize	|
--track_cn			runs IGV_TRACK_CNS in HUMAN_WGS workflow
--track_sv			runs IGV_TRACK_VCF_SV_{jabba|manta}
--track_read			runs IGV_TRACK_READ
"""
    
    exit 0
}

include {
	parse_stub_json
} from "./modules/stub"

// Optionally load json map to control the behaviour of stubs (cp vs touch)
stub_json_map = params.stub_json && ( params.stub_json.toString ().toLowerCase ().endsWith ("js") || params.stub_json.toString ().toLowerCase ().endsWith ("json") ) ? parse_stub_json (params.stub_json) : null

include {
	extract_data;
	file_has_extension
} from "./modules/input"

include {
	PREPARE_GENOME;
	GENOME_ANNOTATION
} from "./modules/local/subworkflow/genome"

include {
	MAP as MAPPER
	REMAP as REMAPPER
} from "./modules/local/subworkflow/remap"

include {
	MANTA
} from "./modules/local/subworkflow/manta" addParams (stub_json_map: stub_json_map)

include {
	STRELKA
} from "./modules/local/subworkflow/strelka"

include {
	MUTECT;
	MUTECT_ANNOTATE;
	MUTECT_RESULT;
	MUTECT_RESULT_RARE
} from "./modules/local/subworkflow/mutect" addParams (stub_json_map: stub_json_map)

include {
	DELLY
} from "./modules/local/subworkflow/delly"

include {
	CNV_KIT;
	CNV_KIT_SELECT_SAMPLE;
	CNV_KIT_TARGET_BED;
	CNV_KIT_COVERAGE;
	CNV_KIT_FIX;
	CNV_KIT_SEGMENT;
	CNV_KIT_PON
} from "./modules/local/subworkflow/cnv-kit" addParams (stub_json_map: stub_json_map)

include {
	HMM_COPY
} from "./modules/local/subworkflow/hmm-copy"

include {
	LOH
} from "./modules/local/subworkflow/loh"

include {
	MSI_SENSOR
} from "./modules/local/subworkflow/msi-sensor"

include {
	BUBBLE_TREE
} from "./modules/local/subworkflow/bubble-tree"

include {
	JABBA
} from "./modules/local/subworkflow/jabba"

include {
	IGV_TRACK_READ;
	IGV_TRACK_CNR as IGV_TRACK_CNR_cnv_kit;
	IGV_TRACK_CNR as IGV_TRACK_CNR_hmm_copy;
	IGV_TRACK_CNS;
	IGV_TRACK_CNS as IGV_TRACK_CNS_cnv_kit;
	IGV_TRACK_VCF_SV as IGV_TRACK_VCF_SV_jabba;
	IGV_TRACK_VCF_SV as IGV_TRACK_VCF_SV_manta
} from "./modules/local/subworkflow/igv-track"

include {
	COHORT_QC as COHORT_QC_human;
	COHORT_QC as COHORT_QC_mouse
} from "./modules/local/subworkflow/qc"

tsv_path = null
pon_bed_path = null
pon_tsv_path = null

ch_input_sample = Channel.empty ()

// pipeline info
log.info """\
MoCaSeq cancer genome sequencing pipeline for mouse and human
=============================================================
Core parameters for run
cache_base		: ${params.cache_base}
script_base		: ${params.script_base}
genomes_base		: ${params.genomes_base}
genome_build.human	: ${params.genome_build.human}
genome_build.mouse	: ${params.genome_build.mouse}
genomes			: ${params.genomes}
genome_annotations	: ${params.genome_annotations}
input			: ${params.input}
output_base		: ${params.output_base}
""".stripIndent()
// dry_init option stops before executing any steps, but prints run configuration
if (params.dry_init != null) { exit 0 }

// check if we have valid --input
if (params.input == null && params.pon_tsv == null && params.qc_dir == null) {
	  exit 1, "[MoCaSeq] error: --input or --pon_tsv or --qc_dir were not supplied! Please check '--help' or documentation under 'running the pipeline' for details"
}

// Read in files properly from TSV file
if (params.input && (file_has_extension (params.input, "tsv"))) tsv_path = params.input
if (params.pon_bed && (file_has_extension (params.pon_bed, "bed"))) pon_bed_path = params.pon_bed
if (params.pon_tsv && (file_has_extension (params.pon_tsv, "tsv"))) pon_tsv_path = params.pon_tsv

if (tsv_path) {

	tsv_file = file (tsv_path)
	if (tsv_file instanceof List) exit 1, "[MoCaSeq] error: can only accept one TSV file per run."
	if (!tsv_file.exists ()) exit 1, "[MoCaSeq] error: input TSV file could not be found. Does the file exist and is it in the right place? You gave the path: ${params.input}"
	ch_input_sample = extract_data (tsv_path)

}
else if (pon_tsv_path)
{
	tsv_file = file (pon_tsv_path)
	if (tsv_file instanceof List) exit 1, "[MoCaSeq] error: can only accept one TSV file per run."
	if (!tsv_file.exists ()) exit 1, "[MoCaSeq] error: pon_tsv TSV file could not be found. Does the file exist and is it in the right place? You gave the path: ${params.pon_tsv}"
}
else if (params.qc_dir)
{
	qc_dir_file = file (params.qc_dir)
	if (qc_dir_file instanceof List) exit 1, "[MoCaSeq] error: can only accept one QC dir per run."
	if ( ! ( qc_dir_file.exists () && qc_dir_file.isDirectory () ) ) exit 1, "[MoCaSeq] error: qc_dir directory could not be found or was not a directory. You gave the path: ${params.qc_dir}"
}
else exit 1, "[MoCaSeq] error: --input or --pon_tsv file(s) or --qc_dir not supplied or improperly defined, see '--help' flag and documentation under 'running the pipeline' for details."

ch_input_branched = ch_input_sample.branch {
	bam: it["normalBAM"] != null || it["tumorBAM"] != null //These are all BAMs
	map: ( it["normalR1"] != null && it["normalR1"].toString().endsWith (".fastq.gz") ) || ( it["tumorR1"] != null && it["tumorR1"].toString().endsWith (".fastq.gz") ) //Path.endsWith tries to match entire final segment
	remap: ( it["normalR1"] != null && it["normalR1"].toString().endsWith (".bam") ) || ( it["tumorR1"] != null && it["tumorR1"].toString().endsWith (".bam") ) //Path.endsWith tries to match entire final segment
}

ch_input_branched_bam_branched = ch_input_branched.bam.branch {
	human_wgs: it["organism"] == "human" && it["seqType"] == "wgs"
	mouse_wex: it["organism"] == "mouse" && it["seqType"] == "wex"
	other: true
}

ch_input_branched_map_branched = ch_input_branched.map.branch {
	human_wgs: it["organism"] == "human" && it["seqType"] == "wgs"
	mouse_wgs: it["organism"] == "mouse" && it["seqType"] == "wgs"
	other: true
}

ch_input_branched_remap_branched = ch_input_branched.remap.branch {
	human_wgs: it["organism"] == "human" && it["seqType"] == "wgs"
	mouse_wgs: it["organism"] == "mouse" && it["seqType"] == "wgs"
	mouse_wex: it["organism"] == "mouse" && it["seqType"] == "wex"
	other: true
}

ch_input_branched_bam_branched.other.view { "[MoCaSeq] error: Failed to find matching workflow (organism & seqType) for input bam:\n${it}" }
ch_input_branched_map_branched.other.view { "[MoCaSeq] error: Failed to find matching workflow (organism & seqType) for input map:\n${it}" }
ch_input_branched_remap_branched.other.view { "[MoCaSeq] error: Failed to find matching workflow (organism & seqType) for input remap:\n${it}" }

workflow HUMAN_WGS
{
	main:
	PREPARE_GENOME (params.genome_build.human)
	GENOME_ANNOTATION (params.genome_build.human)

	ch_bam = ch_input_branched_bam_branched.human_wgs

	MANTA (params.genome_build.human, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.interval_bed, ch_bam, GENOME_ANNOTATION.out.snpeff_version)
	STRELKA (params.genome_build.human, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.interval_bed, ch_bam, MANTA.out.indel)
	MUTECT (params.genome_build.human, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out._chrom_n, ch_bam)
	MUTECT_ANNOTATE (params.genome_build.human, PREPARE_GENOME.out.fasta, MUTECT.out.full, GENOME_ANNOTATION.out.snpeff_version, GENOME_ANNOTATION.out.all_vcf, GENOME_ANNOTATION.out.common_vcf, GENOME_ANNOTATION.out.dbnsfp, GENOME_ANNOTATION.out.sift_sources, GENOME_ANNOTATION.out.sift_fields)
	MUTECT_RESULT_RARE (params.genome_build.human, MUTECT_ANNOTATE.out.post_process, MUTECT_ANNOTATE.out.result, GENOME_ANNOTATION.out.cgc, GENOME_ANNOTATION.out.tru_sight)
	DELLY (params.genome_build.human, PREPARE_GENOME.out.fasta, ch_bam)
	HMM_COPY (params.genome_build.human, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.interval_bed, GENOME_ANNOTATION.out.gc_wig, GENOME_ANNOTATION.out.map_wig, ch_bam)
	LOH (params.genome_build.human, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.fasta_index, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.chrom_names_auto, PREPARE_GENOME.out.interval_bed, GENOME_ANNOTATION.out.alt_haplotype, GENOME_ANNOTATION.out.centromere, GENOME_ANNOTATION.out.mappability, MUTECT.out.result)
	MSI_SENSOR (params.genome_build.human, GENOME_ANNOTATION.out.micro_satellite, ch_bam)

	if ( params.pon_dir == null )
	{
		CNV_KIT (params.genome_build.human, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.fasta_index, GENOME_ANNOTATION.out.ref_flat, PREPARE_GENOME.out.interval_bed, GENOME_ANNOTATION.out.gencode_genes_bed, ch_bam)
		BUBBLE_TREE (params.genome_build.human, PREPARE_GENOME.out.chrom_names_auto, HMM_COPY.out.call, LOH.out.result)
		JABBA (params.genome_build.human, PREPARE_GENOME.out.chrom_names, MANTA.out.vcf, HMM_COPY.out.result, BUBBLE_TREE.out.result)

		if ( params.track_cn )
		{
			IGV_TRACK_CNS (params.genome_build.human, CNV_KIT.out.call)
		}
	}
	else
	{
		ch_bam_tumor = ch_bam.filter { it["type"] == "Tumor" }.map { tuple (it, "Tumor", it["tumorBAM"], it["tumorBAI"] ) }
		ch_target_bed = Channel.of ([
			file ("${params.pon_dir}/${params.genome_build.human}_PON/${params.genome_build.human}.target.bed", glob: false),
			file ("${params.pon_dir}/${params.genome_build.human}_PON/${params.genome_build.human}.resolution.json", glob: false)
		]).map {
			def jsonSlurper = new groovy.json.JsonSlurper ()
			def data = jsonSlurper.parse (it[1])
			tuple ( it[0], data["target_avg_size"], data["wgs_depth"] )
		}.first ()

		ch_centre = params.cnv_kit && params.cnv_kit.centre ? Channel.value (params.cnv_kit.centre.tokenize (",")) : Channel.empty ()

		CNV_KIT_COVERAGE (params.genome_build.human, PREPARE_GENOME.out.fasta, ch_target_bed, ch_bam_tumor)
		CNV_KIT_FIX (params.genome_build.human, Channel.fromPath ("${params.pon_dir}/${params.genome_build.human}_PON/${params.genome_build.human}.reference.cnn").first (), CNV_KIT_COVERAGE.out.result)
		CNV_KIT_SEGMENT (params.genome_build.human, PREPARE_GENOME.out.interval_bed, ch_centre, CNV_KIT_FIX.out.cnr)
		BUBBLE_TREE (params.genome_build.human, PREPARE_GENOME.out.chrom_names_auto, CNV_KIT_SEGMENT.out.call, LOH.out.result)
		JABBA (params.genome_build.human, PREPARE_GENOME.out.chrom_names, MANTA.out.vcf, CNV_KIT_SEGMENT.out.result, BUBBLE_TREE.out.result)

		if ( params.track_cn )
		{
			IGV_TRACK_CNR_cnv_kit (params.genome_build.human,  PREPARE_GENOME.out.interval_bed, CNV_KIT_SEGMENT.out.cnr)
			IGV_TRACK_CNR_hmm_copy (params.genome_build.human, PREPARE_GENOME.out.interval_bed, HMM_COPY.out.cnr)
			IGV_TRACK_CNS_cnv_kit (params.genome_build.human, CNV_KIT_SEGMENT.out.call)
		}
		if ( params.track_sv )
		{
			IGV_TRACK_VCF_SV_jabba (params.genome_build.human, JABBA.out.vcf.map { [ it[0], "JaBbA", it[1] ] })
			IGV_TRACK_VCF_SV_manta (params.genome_build.human, MANTA.out.vcf.map { [ it[0], "Manta", it[1] ] })
		}
	}

	if ( params.track_read )
	{
		IGV_TRACK_READ (params.genome_build.human, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.interval_bed, ch_bam)
	}
}

workflow MOUSE_WEX
{
	main:
	PREPARE_GENOME (params.genome_build.mouse)
	GENOME_ANNOTATION (params.genome_build.mouse)

	ch_bam = ch_input_branched_bam_branched.mouse_wex

	MANTA (params.genome_build.mouse, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.interval_bed, ch_bam, GENOME_ANNOTATION.out.snpeff_version)
	STRELKA (params.genome_build.mouse, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.interval_bed, ch_bam, MANTA.out.indel)
	MUTECT (params.genome_build.mouse, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out._chrom_n, ch_bam)
	MUTECT_ANNOTATE (params.genome_build.mouse, PREPARE_GENOME.out.fasta, MUTECT.out.full, GENOME_ANNOTATION.out.snpeff_version, GENOME_ANNOTATION.out.all_vcf, GENOME_ANNOTATION.out.common_vcf, GENOME_ANNOTATION.out.dbnsfp, GENOME_ANNOTATION.out.sift_sources, GENOME_ANNOTATION.out.sift_fields)
	MUTECT_RESULT (params.genome_build.mouse, MUTECT_ANNOTATE.out.post_process, MUTECT_ANNOTATE.out.result, GENOME_ANNOTATION.out.cgc)
	HMM_COPY (params.genome_build.mouse, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.interval_bed, GENOME_ANNOTATION.out.gc_wig, GENOME_ANNOTATION.out.map_wig, ch_bam)
	LOH (params.genome_build.mouse, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.fasta_index, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.chrom_names_auto, PREPARE_GENOME.out.interval_bed, GENOME_ANNOTATION.out.alt_haplotype, GENOME_ANNOTATION.out.centromere, GENOME_ANNOTATION.out.mappability, MUTECT.out.result)
	MSI_SENSOR (params.genome_build.mouse, GENOME_ANNOTATION.out.micro_satellite, ch_bam)
	BUBBLE_TREE (params.genome_build.mouse, PREPARE_GENOME.out.chrom_names_auto, HMM_COPY.out.call, LOH.out.result)
}

workflow HUMAN_MAP {
	main:
	if (params.debug) { println "[MoCaSeq] debug: entered HUMAN_MAP workflow" }
	PREPARE_GENOME (params.genome_build.human)
	GENOME_ANNOTATION (params.genome_build.human)

	MAPPER (params.genome_build.human, PREPARE_GENOME.out.bwa_index, PREPARE_GENOME.out.fasta, GENOME_ANNOTATION.out.common_vcf, ch_input_branched_map_branched.human_wgs)
	if (params.debug) { println "[MoCaSeq] debug: bwa_index from prep genome out" }
	PREPARE_GENOME.out.bwa_index.view()
	REMAPPER (params.genome_build.human, PREPARE_GENOME.out.bwa_index, PREPARE_GENOME.out.fasta, GENOME_ANNOTATION.out.common_vcf, ch_input_branched_remap_branched.human_wgs)
}

workflow MOUSE_MAP {
	main:
	PREPARE_GENOME (params.genome_build.mouse)
	GENOME_ANNOTATION (params.genome_build.mouse)

	MAPPER (params.genome_build.mouse, PREPARE_GENOME.out.bwa_index, PREPARE_GENOME.out.fasta, GENOME_ANNOTATION.out.common_vcf, ch_input_branched_map_branched.mouse_wgs)
	REMAPPER (params.genome_build.mouse, PREPARE_GENOME.out.bwa_index, PREPARE_GENOME.out.fasta, GENOME_ANNOTATION.out.common_vcf, ch_input_branched_remap_branched.mouse_wgs)
}

workflow HUMAN_PON {
	main:
	PREPARE_GENOME (params.genome_build.human)
	GENOME_ANNOTATION (params.genome_build.human)

	if ( pon_bed_path == null )
	{
		// Generate target regions for CNVKit
		ch_bam_normal = ch_input_branched_bam_branched.human_wgs.filter { it["type"] == "Normal" }.map { tuple (it, "Normal", it["normalBAM"], it["normalBAI"] ) }

		if ( params.pon_sample == null )
		{
			// First we need to pick a median coverage bam
			ch_normal_size_lines = ch_bam_normal.map { [ it[0]["sampleName"], params.genome_build.human, it[0]["normalBAM"], it[0]["normalBAI"], it[0]["normalBAM"].size () ].join ("\t") }
			ch_normal_size_tsv = Channel.of ( ["sample", "genome_build", "bam_path", "bai_path", "bam_size"].join ("\t") )
				.concat (ch_normal_size_lines)
				.collectFile (name: "${params.genome_build.human}.normal_sizes.tsv", newLine: true, sort: false, storeDir: "${params.output_base}/${params.genome_build.human}_PON")

			if ( params.pon_dry )
			{
				CNV_KIT_SELECT_SAMPLE (params.genome_build.human, ch_normal_size_tsv)
				ch_bam_normal_sample = ch_bam_normal.map { [ it[0]["sampleName"], it] }.join (CNV_KIT_SELECT_SAMPLE.out.result).map { it[1] }
				CNV_KIT_TARGET_BED (params.genome_build.human, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.chrom_names, ch_bam_normal_sample)
				CNV_KIT_COVERAGE (params.genome_build.human, PREPARE_GENOME.out.fasta, CNV_KIT_TARGET_BED.out.result, ch_bam_normal)

				ch_normal_coverage_lines = CNV_KIT_COVERAGE.out.result.map { [it[0]["sampleName"], params.genome_build.human, it[2], it[3]].join ("\t") }
				ch_normal_coverage_tsv = Channel.of ( ["sample", "genome_build", "resolution", "normal_cov"].join ("\t")  )
					.concat (ch_normal_coverage_lines)
					.collectFile (name: "${params.genome_build.human}.normal_coverage_file_paths.tsv", newLine: true, sort: false, storeDir: "${params.output_base}/${params.genome_build.human}_PON")

				CNV_KIT_PON (params.genome_build.human, PREPARE_GENOME.out.fasta, ch_normal_coverage_tsv)
			}
		}
		else
		{
			ch_bam_normal_sample = ch_bam_normal.first { it[0]["sampleName"] == params.pon_sample }
			CNV_KIT_TARGET_BED (params.genome_build.human, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.chrom_names, ch_bam_normal_sample)
		}
	}
	else
	{
		if ( pon_tsv_path == null )
		{
			if ( params.pon_resolution == null ) exit 1, "[MoCaSeq] error: You must also supply --pon_resolution when you give --pon_tsv (to ensure consistent names for the coverage .cnn files)"
			ch_bam_normal = ch_input_branched_bam_branched.human_wgs.filter { it["type"] == "Normal" }.map { tuple (it, "Normal", it["normalBAM"], it["normalBAI"] ) }
			ch_target_bed = Channel.value ( [ file (pon_bed_path, glob: false), params.pon_resolution as long, 0 ] )

			CNV_KIT_COVERAGE (params.genome_build.human, PREPARE_GENOME.out.fasta, ch_target_bed, ch_bam_normal)

			ch_normal_coverage_lines = CNV_KIT_COVERAGE.out.result.map { [it[0]["sampleName"], params.genome_build.human, it[2], it[3]].join ("\t") }
			ch_normal_coverage_tsv = Channel.of ( ["sample", "genome_build", "resolution", "normal_cov"].join ("\t")  )
				.concat (ch_normal_coverage_lines)
				.collectFile (name: "${params.genome_build.human}.normal_coverage_file_paths.tsv", newLine: true, sort: false, storeDir: "${params.output_base}/${params.genome_build.human}_PON")

			if ( params.pon_dry )
			{
				CNV_KIT_PON (params.genome_build.human, PREPARE_GENOME.out.fasta, ch_normal_coverage_tsv)
			}
		}
		else
		{
			CNV_KIT_PON (params.genome_build.human, PREPARE_GENOME.out.fasta, Channel.fromPath (pon_tsv_path))
		}
	}
}

workflow MOUSE_PON {
	main:
	PREPARE_GENOME (params.genome_build.mouse)
	GENOME_ANNOTATION (params.genome_build.mouse)

	if ( pon_tsv_path == null )
	{
		ch_bam_normal = ch_input_branched_bam_branched.mouse_wex.filter { it["type"] == "Normal" }.map { tuple (it, "Normal", it["normalBAM"], it["normalBAI"] ) }

		CNV_KIT_COVERAGE (params.genome_build.mouse, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.interval_bed, ch_bam_normal)

		ch_normal_coverage_lines = CNV_KIT_COVERAGE.out.result.map { [it[0]["sampleName"], params.genome_build.mouse, it[2], it[3]].join ("\t") }
		ch_normal_coverage_tsv = Channel.of ( ["sample", "genome_build", "resolution", "normal_cov"].join ("\t")  )
			.concat (ch_normal_coverage_lines)
			.collectFile (name: "${params.genome_build.mouse}.normal_coverage_file_paths.tsv", newLine: true, sort: false, storeDir: "${params.output_base}/${params.genome_build.mouse}_PON")

		if ( params.pon_dry )
		{
			//CNV_KIT_PON (params.genome_build.mouse, PREPARE_GENOME.out.chrom_names, GENOME_ANNOTATION.out.par_interval_bed, ch_normal_coverage_tsv)
		}
	}
	else
	{
		//CNV_KIT_PON (params.genome_build.mouse, PREPARE_GENOME.out.chrom_names, GENOME_ANNOTATION.out.par_interval_bed, Channel.fromPath (pon_tsv_path))
	}
}

workflow {
	if ( params.genome_build.human ) { HUMAN_WGS () }
	if ( params.genome_build.mouse ) { MOUSE_WEX () }
}

// Run using -entry MAP
workflow MAP {
	if (params.debug) { println "[MoCaSeq] debug: entered MAP worfklow" }
	if ( params.genome_build.human ) { HUMAN_MAP () }
	if ( params.genome_build.mouse ) { MOUSE_MAP () }
}

// Run using -entry PON
workflow PON {
	if ( params.genome_build.human ) { HUMAN_PON () }
	if ( params.genome_build.mouse ) { MOUSE_PON () }
}


// Run using -entry QC
workflow QC {
	COHORT_QC_human (params.genome_build.human, Channel.fromPath (params.qc_dir))
	COHORT_QC_mouse (params.genome_build.mouse, Channel.fromPath (params.qc_dir))
}

