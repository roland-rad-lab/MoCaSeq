
include { dry_clean_detergent;dry_clean } from "../software/dry-clean/main"


workflow DRY_CLEAN {

	take:
		ch_interval
		pon_dir_path
		ch_tumor_coverage

	main:
	ch_interval_csv_string = ch_interval.toList ().map { it.join (",") }

	def dir_pon = file (pon_dir_path, glob: false)
	if ( !dir_pon.exists ()) exit 1, "[MoCaSeq] error: PON directory '${dir_pon}' not found."
	if ( !dir_pon.isDirectory ()) exit 1, "[MoCaSeq] error: PON directory required '${dir_pon}' is not a directory."

	ch_pon = Channel.of ([dir_pon.resolve ("normal_table.rds"), dir_pon.resolve ("germline.markers.rds"), dir_pon.resolve ("detergent.rds") ])

	dry_clean (ch_interval_csv_string, ch_pon, ch_tumor_coverage)

	emit:
		tsv = dry_clean.out.result
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

