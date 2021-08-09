
include { dry_clean_detergent;dry_clean } from "../software/dry-clean/main"


workflow DRY_CLEAN {

	take:
		genome_build
		ch_interval
		pon_dir_path
		ch_tumor_coverage

	main:
	ch_interval_csv_string = ch_interval.toList ().map { it.join (",") }

	def dir_pon = file (pon_dir_path, glob: false)
	if ( !dir_pon.exists ()) exit 1, "[MoCaSeq] error: PON directory '${dir_pon}' not found."
	if ( !dir_pon.isDirectory ()) exit 1, "[MoCaSeq] error: PON directory required '${dir_pon}' is not a directory."

	ch_pon = Channel.of ([dir_pon.resolve ("${genome_build}.normal_table.rds"), dir_pon.resolve ("${genome_build}.germline.markers.filtered.rds"), dir_pon.resolve ("${genome_build}.detergent.rds") ]).first ()

	dry_clean (genome_build, ch_interval_csv_string, ch_pon, ch_tumor_coverage)

	emit:
		tsv = dry_clean.out.result
		cnr = dry_clean.out.cnr
}

workflow DRY_CLEAN_PON {

	take:
		genome_build
		ch_interval
		ch_par_interval_bed
		ch_normal_coverage_tsv

	main:
	ch_interval_csv_string = ch_interval.toList ().map { it.join (",") }

	ch_normal_coverage_tsv_filtered = ch_normal_coverage_tsv.splitCsv (header: true, sep: "\t")
		.filter {
			if ( !it.containsKey ("genome_build") ) exit 1, "[MoCaSeq] error: Invalid lines in normal_coverage_tsv, no 'genome_build' column '${it}'"
			return it["genome_build"] == genome_build
		}.tap { ch_paths }
		.map { jt ->
			def f = file (jt["normal_cov"],glob: false)
			// DO NOT MODIFY jt
			def header = ["sample", "genome_build", "normal_cov"]
			def ldata = [jt["sample"], jt["genome_build"], f.name]
			[ header.join ("\t"), ldata.join ("\t"), ""].join ("\n")
		}
		.collectFile (keepHeader: true, skip: 1)

	dry_clean_detergent (genome_build, ch_interval_csv_string, ch_par_interval_bed, ch_paths.map { file (it["normal_cov"], glob: false) }.collect (), ch_normal_coverage_tsv_filtered)
}

