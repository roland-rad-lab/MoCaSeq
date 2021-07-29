
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

	ch_pon = Channel.of ([dir_pon.resolve ("${genome_build}.normal_table.rds"), dir_pon.resolve ("${genome_build}.germline.markers.rds"), dir_pon.resolve ("${genome_build}.detergent.rds") ])

	dry_clean (ch_interval_csv_string, ch_pon, ch_tumor_coverage)

	emit:
		tsv = dry_clean.out.result
}

workflow DRY_CLEAN_PON {

	take:
		genome_build
		ch_interval
		ch_par_interval_bed
		ch_normal_coverage_tsv

	main:
	ch_interval_csv_string = ch_interval.toList ().map { it.join (",") }

	ch_normal_coverage_tsv_filtered = ch_normal_coverage_tsv.splitText (keepHeader: true)
		.filter {
			def lines = it.tokenize ("\n")
			if ( lines.size () < 2 ) exit 1, "[MoCaSeq] error: Invalid lines in normal_coverage_tsv '${it}'"
			def m = [lines[0],lines[1]].transpose().collectEntries()
			if ( !m.containsKey ("genome_build") ) exit 1, "[MoCaSeq] error: Invalid lines in normal_coverage_tsv, no 'genome_build' column '${it}'"
			return m["genome_build"] == genome_build
		}
		.tap { ch_have_work }
		.collectFile (keepHeader: true, skip: 1)

	dry_clean_detergent (genome_build, ch_interval_csv_string, ch_par_interval_bed, ch_have_work.count ().filter { it > 0 }, ch_normal_coverage_tsv)
}

