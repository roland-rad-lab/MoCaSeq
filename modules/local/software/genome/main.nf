
boolean valid_uri (String maybe_url)
{
	try {
		java.net.URI.create (maybe_url);
		true
	} catch ( IllegalArgumentException | NullPointerException e ) {
		false
	}
}

void download_uri (Path path, java.net.URI uri)
{
	try {
		java.nio.file.Files.createFile (path)
		println ("download '${uri}' to '${path}'")
		path.toFile ().withOutputStream { out -> uri.toURL ().withInputStream { in -> out << in; } }
	} catch ( java.nio.file.FileAlreadyExistsException e) {
		println ("skipping '${path}', already exists")
	} catch ( IOException e) {
		exit 1, "Failed to download '${uri}'"
	}
}

process cache_genome_url {

	input:
		val (genome_build)
		val (reference)
		val (extension_list)

	output:
		val (path_cached), emit: result

	exec:
		def reference_local = reference;
		if ( extension_list.size () == 0 ) { exit 1, "[MoCaSeq] error: At least one file extension is required to append to '${reference_local}'" }
		if ( valid_uri (reference_local) )
		{
			def reference_uri = java.net.URI.create (reference_local)
			if ( reference_uri.scheme && reference_uri.scheme == "http" ) { exit 1, "[MoCaSeq] error: No support for http downloads use https, url requested was '${reference_local}'" }
			if ( reference_uri.scheme && reference_uri.scheme == "https" && reference_uri.path )
			{
				def cache_dir = workDir.resolve ("${genome_build}_cache");
				java.nio.file.Files.createDirectories (cache_dir);
				def file_name_base = reference_uri.path.tokenize ('/')[-1]
				path_cached = cache_dir.resolve ( [file_name_base, extension_list[0]].findAll { it }.join (".") )

				extension_list.findAll { it != null }.each {
					def file_name_extended = [file_name_base, it].findAll { jt -> jt != "" }.join (".")
					if ( !valid_uri ( file_name_extended ) ) { exit 1, "[MoCaSeq] error: The extension '${file_name_extended}' was not a valid URI" }
					download_uri (cache_dir.resolve (file_name_extended), reference_uri.resolve ( java.net.URI.create (file_name_extended)))
				}
			}
			else
			{
				path_cached = reference_local.resolveSibling ( [reference_local.fileName.toString (), extension_list[0]].join (".") )
			}
		}
		else
		{
			path_cached = reference_local.resolveSibling ( [reference_local.fileName.toString (), extension_list[0]].join (".") )
		}
}

process interval_bed {

	input:
		val (dict)
		path (interval_list)

	output:
		tuple path ("intervals.bed.gz"), path ("intervals.bed.gz.tbi"), emit: result

	script:
	"""#!/usr/bin/env Rscript

library (dplyr)
library (stringr)

dict_file_path = "${dict}"
intervals_file_path = "interval_names.tsv"
output_bed_file_path = "intervals.bed.gz"

data_dict <- read.table (file=dict_file_path,sep="\\t",stringsAsFactors=F,header=F,skip=1)
names (data_dict) <- c("line_type", "sequence_name_raw", "sequence_length_raw", "md5", "url")
head (data_dict)

data_seq_lengths <- data_dict %>%
  mutate (sequence_name=stringr::str_split_fixed (sequence_name_raw,":",2)[,2]) %>%
  mutate (sequence_length=as.numeric (stringr::str_split_fixed (sequence_length_raw,":",2)[,2])) %>%
  select (sequence_name,sequence_length) %>%
  data.frame

head (data_seq_lengths)

data_intervals <- read.table (file=intervals_file_path,sep="\\t",stringsAsFactors=F,header=F)
names (data_intervals) <- c("sequence_name")

data_output <- data_intervals %>%
	inner_join (data_seq_lengths,by="sequence_name") %>%
	mutate (start=0) %>%
	select (sequence_name,start,sequence_length) %>%
	data.frame

output_bed_file <-pipe (paste ("bgzip -c >",output_bed_file_path), "w")
write.table (data_output,file=output_bed_file,sep="\\t",row.names=F,col.names=F,quote=F)
close (output_bed_file)
system2 ("tabix",args=c("-p", "vcf", output_bed_file_path),wait=T)
"""

}

process interval_bed_intersect {

	input:
		val (bed_a)
		path (bed_b)
		val (flags)

	output:
		tuple path ("intervals.intersection.bed.gz"), path ("intervals.intersection.bed.gz.tbi"), emit: result

	script:
	"""#!/usr/bin/env bash
source ${params.script_base}/file_handling.sh
temp_file_b=\$(moc_mktemp_file .)
trap "rm \${temp_file_b}" EXIT

extract_if_zip ${bed_b} bed_b_extracted \${temp_file_b}

bedtools intersect -a ${bed_a} -b \${bed_b_extracted} ${flags} | bgzip -c > intervals.intersection.bed.gz
tabix -p bed intervals.intersection.bed.gz

	"""
}


