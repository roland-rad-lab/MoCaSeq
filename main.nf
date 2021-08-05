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
	CNV_KIT;
	CNV_KIT_SEGMENT
} from "./modules/local/subworkflow/cnv-kit"

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
	IGV_TRACK_CNS;
	IGV_TRACK_CNS as IGV_TRACK_CNS_dryclean;
	IGV_TRACK_RDS
} from "./modules/local/subworkflow/igv-track"

include {
	FRAG_COUNTER
} from "./modules/local/subworkflow/frag-counter"

include {
	DRY_CLEAN;
	DRY_CLEAN_PON
} from "./modules/local/subworkflow/dry-clean"

tsv_path = null
pon_tsv_path = null


ch_input_sample = Channel.empty ()


// check if we have valid --input
if (params.input == null && params.pon_tsv == null) {
	  exit 1, "[MoCaSeq] error: --input or --pon_tsv was not supplied! Please check '--help' or documentation under 'running the pipeline' for details"
}

// Read in files properly from TSV file
if (params.input && (file_has_extension (params.input, "tsv"))) tsv_path = params.input
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
} else exit 1, "[MoCaSeq] error: --input or --pon_tsv file(s) not correctly not supplied or improperly defined, see '--help' flag and documentation under 'running the pipeline' for details."

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
	MSI_SENSOR (GENOME_ANNOTATION.out.micro_satellite, ch_bam)

	if ( params.pon_dir == null )
	{
		BUBBLE_TREE (PREPARE_GENOME.out.chrom_names_auto, HMM_COPY.out.tsv, LOH.out.result)
		JABBA (MANTA.out.basic, HMM_COPY.out.tsv, BUBBLE_TREE.out.result)
	}
	else
	{
		FRAG_COUNTER (params.genome_build.human, PREPARE_GENOME.out.chrom_names, GENOME_ANNOTATION.out.gc_wig, GENOME_ANNOTATION.out.map_wig, ch_bam)
		DRY_CLEAN (params.genome_build.human, PREPARE_GENOME.out.chrom_names, params.pon_dir, FRAG_COUNTER.out.result)
		CNV_KIT_SEGMENT ("dryclean", DRY_CLEAN.out.cnr)
		BUBBLE_TREE (PREPARE_GENOME.out.chrom_names_auto, CNV_KIT_SEGMENT.out.tsv, LOH.out.result)
		JABBA (MANTA.out.basic, CNV_KIT_SEGMENT.out.tsv, BUBBLE_TREE.out.result)

		if ( params.track_cn )
		{
			IGV_TRACK_RDS (PREPARE_GENOME.out.interval_bed, "fragCounter", FRAG_COUNTER.out.result)
			IGV_TRACK_CNS_dryclean ("dryclean-CNVKit", CNV_KIT_SEGMENT.out.cns)
		}
	}

	if ( params.track_read )
	{
		IGV_TRACK_READ (PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.interval_bed, ch_bam)
	}
	if ( params.track_cn )
	{
		IGV_TRACK_CNS ("CNVKit", CNV_KIT.out.cns)
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
	MUTECT (PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out._chrom_n, ch_bam)
	HMM_COPY (PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.interval_bed, GENOME_ANNOTATION.out.gc_wig, GENOME_ANNOTATION.out.map_wig, ch_bam)
	LOH (PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.fasta_index, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.interval_bed, MUTECT.out.result)
	MSI_SENSOR (GENOME_ANNOTATION.out.micro_satellite, ch_bam)
	BUBBLE_TREE (PREPARE_GENOME.out.chrom_names_auto, HMM_COPY.out.tsv, LOH.out.result)
}

workflow HUMAN_PON {
	main:
	PREPARE_GENOME (params.genome_build.human)
	GENOME_ANNOTATION (params.genome_build.human)

	if ( pon_tsv_path == null )
	{
		ch_bam_normal = ch_input_branched_bam_branched.human_wgs.map { tuple (it, "Normal", it["normalBAM"], it["normalBAI"] ) }

		FRAG_COUNTER (params.genome_build.human, PREPARE_GENOME.out.chrom_names, GENOME_ANNOTATION.out.gc_wig, GENOME_ANNOTATION.out.map_wig, ch_bam_normal)

		ch_normal_coverage_lines = FRAG_COUNTER.out.result.map { [it[0]["sampleName"], params.genome_build.human, it[2]].join ("\t") }
		ch_normal_coverage_tsv = Channel.of ( ["sample", "genome_build", "normal_cov"].join ("\t")  )
			.concat (ch_normal_coverage_lines)
			.collectFile (name: "${params.genome_build.human}.normal_coverage_file_paths.tsv", newLine: true, sort: false, storeDir: "${params.output_base}/PON")

		if ( params.pon_dry )
		{
			DRY_CLEAN_PON (params.genome_build.human, PREPARE_GENOME.out.chrom_names, GENOME_ANNOTATION.out.par_interval_bed, ch_normal_coverage_tsv)
		}
	}
	else
	{
		DRY_CLEAN_PON (params.genome_build.human, PREPARE_GENOME.out.chrom_names, GENOME_ANNOTATION.out.par_interval_bed, Channel.fromPath (pon_tsv_path))
	}
}

workflow MOUSE_PON {
	main:
	PREPARE_GENOME (params.genome_build.mouse)
	GENOME_ANNOTATION (params.genome_build.mouse)

	if ( pon_tsv_path == null )
	{
		ch_bam_normal = ch_input_branched_bam_branched.mouse_wex.map { tuple (it, "Normal", it["normalBAM"], it["normalBAI"] ) }

		FRAG_COUNTER (params.genome_build.mouse, PREPARE_GENOME.out.chrom_names, GENOME_ANNOTATION.out.gc_wig, GENOME_ANNOTATION.out.map_wig, ch_bam_normal)

		ch_normal_coverage_lines = FRAG_COUNTER.out.result.map { [it[0]["sampleName"], params.genome_build.mouse, it[2]].join ("\t") }
		ch_normal_coverage_tsv = Channel.of ( ["sample", "genome_build", "normal_cov"].join ("\t")  )
			.concat (ch_normal_coverage_lines)
			.collectFile (name: "${params.genome_build.mouse}.normal_coverage_file_paths.tsv", newLine: true, sort: false, storeDir: "${params.output_base}/PON")

		if ( params.pon_dry )
		{
			DRY_CLEAN_PON (params.genome_build.mouse, PREPARE_GENOME.out.chrom_names, GENOME_ANNOTATION.out.par_interval_bed, ch_normal_coverage_tsv)
		}
	}
	else
	{
		DRY_CLEAN_PON (params.genome_build.mouse, PREPARE_GENOME.out.chrom_names, GENOME_ANNOTATION.out.par_interval_bed, Channel.fromPath (pon_tsv_path))
	}
}

workflow {
	HUMAN_WGS ()
	MOUSE_WEX ()
}

// Run using -entry PON
workflow PON {
	HUMAN_PON ()
	MOUSE_PON ()
}


