#!/usr/bin/env nextflow

if (params.help){
	println """\
MoCaSeq test wrapper for HUMAN_WGS:HMM_COPY:hmm_copy_wig_tumor
======================================================================================
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
if (params.dry_init) { exit 0 }

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
	if (params.debug) {
		println "[MoCaSeq] debug: ch_input_sample:"
		ch_input_sample.view()
	}
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

ch_input_branched_bam_branched.other.view { "[MoCaSeq] error: Failed to find matching workflow (organism & seqType) for input bam:\n${it}" }

workflow HUMAN_WGS
{
	main:
	if (params.debug) { println "[MoCaSeq] debug: entered HUMAN_WGS worfklow" }
	
	PREPARE_GENOME (params.genome_build.human)
	GENOME_ANNOTATION (params.genome_build.human)
	
	ch_bam = ch_input_branched_bam_branched.human_wgs

	// debug the ch_bam, to see if the data is parsed correctly
	if (params.debug) {
		println "[MoCaSeq] debug: pre HMM_COPY chanel values"
		println "${params.genome_build.human}"
		PREPARE_GENOME.out.chrom_names.view()
		PREPARE_GENOME.out.interval_bed.view()
		GENOME_ANNOTATION.out.gc_wig.view()
		GENOME_ANNOTATION.out.map_wig.view()
		ch_bam.view()
	}

	HMM_COPY (params.genome_build.human, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.interval_bed, GENOME_ANNOTATION.out.gc_wig, GENOME_ANNOTATION.out.map_wig, ch_bam)
}

workflow {
	if (params.debug) { println "[MoCaSeq] debug: entered main worfklow" }
	if ( params.genome_build.human ) { HUMAN_WGS () }
	if ( params.genome_build.mouse ) { MOUSE_WEX () }
}

