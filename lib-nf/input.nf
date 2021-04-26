#!/usr/bin/env nextflow

// Check file extension
def file_has_extension (it, extension)
{
	it.toString ().toLowerCase ().endsWith (extension.toLowerCase ())
}

// Check if a row has the expected number of columns
def row_check_column_n (row, number)
{
	if (row.size () != number) exit 1, "[MoCaSeq] error:  Invalid TSV input - malformed row (e.g. missing column) in ${row}, see '--help' flag and documentation under 'running the pipeline' for more information"
	return true
}

// Channelling the TSV file containing FASTQ or BAM 
def extract_data (tsv_file)
{
	Channel.fromPath (tsv_file)
		.splitCsv (header: true, sep: '\t')
		.dump (tag:'tsv_extract')
		.map { row ->
			def expected_keys = ['Sample_Name', 'Library_ID', 'Lane', 'Colour_Chemistry', 'SeqType', 'Organism', 'Strandedness', 'UDG_Treatment', 'R1', 'R2', 'BAM']
			if ( !row.keySet ().containsAll (expected_keys) ) exit 1, "[MoCaSeq] error: Invalid TSV input - malformed column names. Please check input TSV. Column names should be: ${expected_keys.join(", ")}"

			row_check_column_n (row, 11)
		}
}

