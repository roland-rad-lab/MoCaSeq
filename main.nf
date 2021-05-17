#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
	extract_data;
	file_has_extension
} from "./modules/input"

include {
	PREPARE_GENOME;
	GENOME_ANNOTATION
} from "./modules/local/subworkflow/genome"

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

ch_input_branched = ch_input_sample.branch {
	bam: it["normalBAM"] != 'NA' //These are all BAMs
}

//Removing R1/R2 in case of BAM input
ch_input_branched_bam_branched = ch_input_branched.bam.branch {
	human: it["organism"] == "human"
	other: true
}

ch_input_branched_bam_branched.other.view { "[MoCaSeq] error: Failed to find matching workflow for input bam:\n${it}" }

workflow
{
	main:
	PREPARE_GENOME (params.genome_build.human)
	GENOME_ANNOTATION (params.genome_build.human)
	MANTA (PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.interval_bed, ch_input_branched_bam_branched.human)
	STRELKA (PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.interval_bed, ch_input_branched_bam_branched.human, MANTA.out.indel)
	MUTECT (PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out._chrom_n, ch_input_branched_bam_branched.human)
	DELLY (PREPARE_GENOME.out.fasta, ch_input_branched_bam_branched.human)
	CNV_KIT (PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.fasta_index_flat, PREPARE_GENOME.out.interval_bed, GENOME_ANNOTATION.out.gencode_genes_bed, ch_input_branched_bam_branched.human)
	HMM_COPY (PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.interval_bed, GENOME_ANNOTATION.out.gc_wig, GENOME_ANNOTATION.out.map_wig, ch_input_branched_bam_branched.human)
	STRUCTURAL_VARIATION (PREPARE_GENOME.out.interval_bed, MANTA.out.result, DELLY.out.result)
}


