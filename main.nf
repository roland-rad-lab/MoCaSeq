#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
	extract_data;
	file_has_extension
} from "./modules/input"

include {
	parse_stub_json
} from "./modules/stub"

include {
	PREPARE_GENOME;
	GENOME_ANNOTATION
} from "./modules/local/subworkflow/genome"

include {
	REMAP
} from "./modules/local/subworkflow/remap"

include {
	MANTA
} from "./modules/local/subworkflow/manta"

include {
	STRELKA
} from "./modules/local/subworkflow/strelka"

include {
	MUTECT
} from "./modules/local/subworkflow/mutect"

include {
	DELLY
} from "./modules/local/subworkflow/delly"

include {
	CNV_KIT
} from "./modules/local/subworkflow/cnv-kit"

include {
	HMM_COPY
} from "./modules/local/subworkflow/hmm-copy"

include {
	LOH
} from "./modules/local/subworkflow/loh"

include {
	BUBBLE_TREE
} from "./modules/local/subworkflow/bubble-tree"

include {
	JABBA
} from "./modules/local/subworkflow/jabba"

include {
	IGV_TRACK_READ;
	IGV_TRACK_CN
} from "./modules/local/subworkflow/igv-track"

tsv_path = null


ch_input_sample = Channel.empty ()


// check if we have valid --input
if (params.input == null) {
	  exit 1, "[MoCaSeq] error: --input was not supplied! Please check '--help' or documentation under 'running the pipeline' for details"
}

// Read in files properly from TSV file
if (params.input && (file_has_extension (params.input, "tsv"))) tsv_path = params.input


if (tsv_path) {

	tsv_file = file (tsv_path)
	if (tsv_file instanceof List) exit 1, "[MoCaSeq] error: can only accept one TSV file per run."
	if (!tsv_file.exists ()) exit 1, "[MoCaSeq] error: input TSV file could not be found. Does the file exist and is it in the right place? You gave the path: ${params.input}"
	ch_input_sample = extract_data (tsv_path)

} else exit 1, "[MoCaSeq] error: --input file(s) not correctly not supplied or improperly defined, see '--help' flag and documentation under 'running the pipeline' for details."

// Optionally load json map to control the behaviour of stubs (cp vs touch)
if (params.stub_json && ( file_has_extension (params.stub_json, "js") || file_has_extension (params.stub_json, "json") ) ) {
	params.stub_json_map = parse_stub_json (params.stub_json)
}

ch_input_branched = ch_input_sample.branch {
	bam: it["normalBAM"] != null //These are all BAMs
	remap: it["normalR1"] != null && it["normalR1"].toString().endsWith (".bam") //Path.endsWith tries to match entire final segment
}

ch_input_branched_bam_branched = ch_input_branched.bam.branch {
	human_wgs: it["organism"] == "human" && it["seqType"] == "wgs"
	mouse_wex: it["organism"] == "mouse" && it["seqType"] == "wex"
	other: true
}

ch_input_branched_remap_branched = ch_input_branched.remap.branch {
	human_wgs: it["organism"] == "human" && it["seqType"] == "wgs"
	mouse_wex: it["organism"] == "mouse" && it["seqType"] == "wex"
	other: true
}

ch_input_branched_bam_branched.other.view { "[MoCaSeq] error: Failed to find matching workflow (organism & seqType) for input bam:\n${it}" }
ch_input_branched_remap_branched.other.view { "[MoCaSeq] error: Failed to find matching workflow (organism & seqType) for input remap:\n${it}" }

workflow HUMAN_WGS
{
	main:
	PREPARE_GENOME (params.genome_build.human)
	GENOME_ANNOTATION (params.genome_build.human)

	REMAP (PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.dir, GENOME_ANNOTATION.out.common_vcf, ch_input_branched_remap_branched.human_wgs)
	ch_bam = REMAP.out.result.mix (ch_input_branched_bam_branched.human_wgs)

	MANTA (PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.interval_bed, ch_bam)
	STRELKA (PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.interval_bed, ch_bam, MANTA.out.indel)
	MUTECT (PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out._chrom_n, ch_bam)
	DELLY (PREPARE_GENOME.out.fasta, ch_bam)
	CNV_KIT (PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.fasta_index_flat, PREPARE_GENOME.out.interval_bed, GENOME_ANNOTATION.out.gencode_genes_bed, ch_bam)
	HMM_COPY (PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.interval_bed, GENOME_ANNOTATION.out.gc_wig, GENOME_ANNOTATION.out.map_wig, ch_bam)
	LOH (PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.fasta_index, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.interval_bed, MUTECT.out.result)
	BUBBLE_TREE (HMM_COPY.out.tsv, LOH.out.result)
	JABBA (MANTA.out.basic, HMM_COPY.out.tsv, BUBBLE_TREE.out.result)

	if ( params.track_read )
	{
		IGV_TRACK_READ (PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.interval_bed, ch_bam)
	}
	if ( params.track_cn )
	{
		IGV_TRACK_CN (CNV_KIT.out.cns_normal, CNV_KIT.out.cns_tumor)
	}
}

workflow MOUSE_WEX
{
	main:
	PREPARE_GENOME (params.genome_build.mouse)
	GENOME_ANNOTATION (params.genome_build.mouse)

	REMAP (PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.dir, GENOME_ANNOTATION.out.common_vcf, ch_input_branched_remap_branched.mouse_wex)
	ch_bam = REMAP.out.result.mix (ch_input_branched_bam_branched.mouse_wex)

	MANTA (PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.interval_bed, ch_bam)
	STRELKA (PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.interval_bed, ch_bam, MANTA.out.indel)
}

workflow {
	HUMAN_WGS ()
	MOUSE_WEX ()
}

