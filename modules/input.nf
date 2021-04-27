#!/usr/bin/env nextflow

// Check file extension
def file_has_extension (it, extension)
{
	it.toString ().toLowerCase ().endsWith (extension.toLowerCase ())
}

// Return file if it exists
def file_from_path (it)
{
	if (!file(it).exists()) exit 1, "[MoCaSeq] error: Cannot find supplied FASTQ or BAM input file. If using input method TSV set to NA if no file required. See '--help' flag and documentation under 'running the pipeline' for more information. Check file: ${it}" 
	return file(it)
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
			def expected_keys = ['Sample_Name', 'Library_ID', 'Lane', 'Colour_Chemistry', 'SeqType', 'Organism', 'Type', 'R1', 'R2', 'BAM']
			if ( !row.keySet ().containsAll (expected_keys) ) exit 1, "[MoCaSeq] error: Invalid TSV input - malformed column names. Please check input TSV. Column names should be: ${expected_keys.join(", ")}"

			row_check_column_n (row, 9)

			if ( row.Sample_Name.isEmpty() ) exit 1, "[MoCaSeq] error: the Sample_Name column is empty. Ensure all cells are filled or contain 'NA' for optional fields. Check row:\n ${row}"
			if ( row.Type.isEmpty () ) exit 1, "[MoCaSeq] error: the Type column is empty. Ensure all cells are filled or contain 'NA' for optional fields. Check row:\n ${row}"
			if ( row.BAM.isEmpty() ) exit 1, "[MoCaSeq] error: the BAM column is empty. Ensure all cells are filled or contain 'NA' for optional fields. Check row:\n ${row}"


			def samplename = row.Sample_Name
			def organism = row.Organism
			def type = row.Type
			def bam = row.BAM.matches('NA') ? 'NA' : file_from_path (row.BAM)

			[ samplename, libraryid, lane, colour, seqtype, organism, type, r1, r2, bam ]
		}.reduce ( [:] ) { accumulator, item ->
			if ( accumulator.containsKey (item[0]) )
			{
				accumulator[item[0]].add (item)
			}
			else
			{
				accumulator[item[0]] = [item]
			}
			accumulator
		}.map { it ->
			it.value.inject ([:]) { accumulator, item ->
				if ( accumulator.size () == 0 )
				{
					accumulator["sampleName"] = it.key
				}
				if ( item[6] == "Normal")
				{
					accumulator["NormalBAM"] = item[8]
				}
				if ( item[6] == "Tumor" )
				{
					accumulator["TumorBAM"] = item[8]
				}
				m
			}
		}
}

