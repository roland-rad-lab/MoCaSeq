#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
	extract_data;
	file_has_extension
} from "./lib-nf/input"

// check if we have valid --reads or --input
if (params.input == null) {
  exit 1, "[nf-core/eager] error: --input was not supplied! Please check '--help' or documentation under 'running the pipeline' for details"
}

// Read in files properly from TSV file
tsv_path = null
if (params.input && (file_has_extensionv(params.input, "tsv"))) tsv_path = params.input

ch_input_sample = Channel.empty()

if (tsv_path) {

	tsv_file = file (tsv_path)
	if (tsv_file instanceof List) exit 1, "[nf-core/eager] error: can only accept one TSV file per run."
	if (!tsv_file.exists ()) exit 1, "[nf-core/eager] error: input TSV file could not be found. Does the file exist and is it in the right place? You gave the path: ${params.input}"
	ch_input_sample = extract_data (tsv_path)

} else exit 1, "[nf-core/eager] error: --input file(s) not correctly not supplied or improperly defined, see '--help' flag and documentation under 'running the pipeline' for details."


