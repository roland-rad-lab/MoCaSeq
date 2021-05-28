

include { sam_to_fastq_paired;fastqc_paired } from "../software/remap/main"

workflow REMAP
{
	take:
		ch_fasta
		ch_data

	main:
		ch_data_branched = ch_data.branch {
			normal_paired: it["normalR1"].endsWith ("bam") && it["normalR1"] == it["normalR2"]
			tumor_paired: it["tumorR1"].endsWith ("bam") && it["tumorR1"] == it["tumorR2"]
			other
		}

		ch_data_branched_normal_paired = ch_data_branched.normal_paired.map { tuple (it, "Normal", it["normalR1"], it["normalR2"] ) }
		ch_data_branched_tumor_paired = ch_data_branched.tumor_paired.map { tuple (it, "Tumor", it["tumorR1"], it["tumorR2"] ) }

		ch_data_expanded_paired = ch_data_branched_normal_paired.mix (ch_data_branched_tumor_paired)

		sam_to_fastq_paired (ch_data_expanded_paired)
		fastqc_paired (sam_to_fastq_paired.out.result)
		foo = Channel.empty ()

	emit:
		result = foo
}

