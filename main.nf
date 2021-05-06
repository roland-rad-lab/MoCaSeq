#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
	extract_data;
	file_has_extension
} from "./modules/input"

include {
	PREPARE_GENOME
} from "./modules/local/subworkflow/genome"

include {
	MUTECT
} from "./modules/local/subworkflow/mutect"


tsv_path = null


ch_input_sample = Channel.empty ()


// check if we have valid --reads or --input
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

ch_branched_input = ch_input_sample.branch {
	bam: it["normalBAM"] != 'NA' //These are all BAMs
}

//Removing R1/R2 in case of BAM input
ch_branched_input_bam_branched = ch_branched_input.bam.branch {
	human: it["organism"] == "Human"
}

workflow
{
	main:
	PREPARE_GENOME (params.genome_build.human)
	MUTECT (PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.chrom_names, ch_branched_input_bam_branched.human)

	ch_delly = MUTECT.out.result.map { it -> tuple ( it[0], it[0]["normalBAM"], it[0]["tumorBAM"] ) }

	DELLY (PREPARE_GENOME.out.fasta, ch_delly)
}


