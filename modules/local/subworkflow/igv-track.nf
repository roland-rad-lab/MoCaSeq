#!/usr/bin/env nextflow

include {
	igv_track_depth as igv_track_depth_normal;
	igv_track_depth as igv_track_depth_tumor;
	igv_track_cnr;
	igv_track_cns;
	igv_track_rds;
	igv_track_vcf_sv
} from "../software/igv-track/main"

workflow IGV_TRACK_READ {

	take:
		genome_build
		ch_interval
		ch_interval_bed
		ch_data

	main:
		ch_interval_space_string = ch_interval.map { it.join (" ") }
		ch_data_normal = ch_data.map { tuple (it, "Normal", it["normalBAM"], it["normalBAI"] ) }
		ch_data_tumor = ch_data.map { tuple (it, "Tumor", it["tumorBAM"], it["tumorBAI"] ) }

		igv_track_depth_normal (genome_build, ch_interval_space_string, ch_interval_bed, ch_data_normal)
		igv_track_depth_tumor (genome_build, ch_interval_space_string, ch_interval_bed, ch_data_tumor)
}

workflow IGV_TRACK_RDS {

	take:
		genome_build
		ch_interval_bed
		coverage_source
		ch_rds

	main:
		igv_track_rds (genome_build, ch_interval_bed, coverage_source, ch_rds)
}

workflow IGV_TRACK_CNR {

	take:
		genome_build
		ch_interval_bed
		ch_cnr

	main:
		igv_track_cnr (genome_build, ch_interval_bed, ch_cnr)
}

workflow IGV_TRACK_CNS {

	take:
		genome_build
		ch_cns

	main:
		igv_track_cns (genome_build, ch_cns)
}

workflow IGV_TRACK_VCF_SV {

	take:
		genome_build
		ch_vcf

	main:
		igv_track_vcf_sv (genome_build, ch_vcf)
}

