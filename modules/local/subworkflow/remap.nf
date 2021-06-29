

include {
	sam_to_fastq_paired;
	fastqc_paired as fastqc_paired_extracted;
	fastqc_paired as fastqc_paired_trimmed;
	trim_paired;
	bwa_mem_paired;
	mark_duplicates_recalibrate
} from "../software/remap/main"

workflow REMAP
{
	take:
		ch_fasta
		ch_dir
		ch_common_vcf
		ch_data

	main:
		ch_data.multiMap { it ->
			normal: tuple (it, "Normal", it["normalR1"], it["normalR2"])
			tumor: tuple (it, "Tumor", it["tumorR1"], it["tumorR2"])
		}.set { ch_data_multi }

		ch_data_multi_branched = ch_data_multi.normal.mix (ch_data_multi.tumor).branch {
			paired: it[2].toString ().endsWith (".bam") && it[2] == it[3]
			other: true
		}

		ch_data_multi_branched.other.view { "[MoCaSeq] error: Failed to find matching REMAP workflow path for input:\n${it}" }

		sam_to_fastq_paired (ch_data_multi_branched.paired.map { tuple (it[0], it[1], it[2] ) })
		fastqc_paired_extracted (sam_to_fastq_paired.out.result)
		trim_paired (fastqc_paired_extracted.out.result)

		trim_paired.out.result.set { ch_trim }

		fastqc_paired_trimmed (ch_trim)
		bwa_mem_paired (ch_fasta, ch_trim)
		mark_duplicates_recalibrate (ch_fasta, ch_common_vcf, bwa_mem_paired.out.result)
		sample = mark_duplicates_recalibrate.out.result.map { [it[0]["sampleName"], it] }
			.groupTuple (size: 2)
			.map { it[1] }
			.map { it ->
				def bam_info = it.inject ([:]) { accumulator, item ->
					accumulator[[item[1].toLowerCase (),"BAM"].join ("")] = item[2]
					accumulator[[item[1].toLowerCase (),"BAI"].join ("")] = item[3]
					accumulator
				}
				def result = it[0][0].clone ()
				result.putAll (bam_info)
				result
			}

	emit:
		result = sample
}

