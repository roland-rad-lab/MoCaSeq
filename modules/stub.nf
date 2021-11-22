#!/usr/bin/env nextflow

def parse_stub_json (json_file_path)
{
	if (json_file_path == null || json_file_path.isEmpty () ) exit 1, "[MoCaSeq] error: No path supplied to parse_stub_json"
	def f = file (json_file_path, glob: false)
	if (!f.exists ()) exit 1, "[MoCaSeq] error: Failed to parse stub json from '${json_file_path}' file does not exist"

	def jsonSlurper = new groovy.json.JsonSlurper ()
	jsonSlurper.parse (f)
}


