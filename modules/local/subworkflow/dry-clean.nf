
include { dry_clean_detergent;dry_clean } from "../software/dry-clean/main"


workflow DRY_CLEAN {

	take:
		pon_dir_path
		ch_tumor_coverage

	main:
	def dir_pon = file (pon_dir_path, glob: false)
	if ( !dir_pon.exists ()) exit 1, "[MoCaSeq] error: PON directory '${dir_pon}' not found."
	if ( !dir_pon.isDirectory ()) exit 1, "[MoCaSeq] error: PON directory required '${dir_pon}' is not a directory."

	def file_normal_table = Paths.get (dir_pon, "normal_table.rds")
	def file_germline_markers = Paths.get (dir_pon, "germline.markers.rds")

	ch_normal_table = Channel.of (file_normal_table)
	ch_file_germline_markers = Channel.of (file_germline_markers)

	dry_clean (ch_normal_table, ch_file_germline_markers, ch_tumor_coverage)
}

workflow DRY_CLEAN_PON {

	take:
		ch_interval
		ch_par_interval_bed
		ch_have_work
		ch_normal_coverage_tsv

	main:
	ch_interval_csv_string = ch_interval.toList ().map { it.join (",") }
	dry_clean_detergent (ch_interval_csv_string, ch_par_interval_bed, ch_have_work, ch_normal_coverage_tsv)
}

