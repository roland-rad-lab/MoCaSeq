

include {
	sam_to_fastq_paired;
	fastqc_paired as fastqc_paired_extracted;
	fastqc_paired as fastqc_paired_trimmed;
	trim_paired;
	bwa_mem_paired;
	mark_duplicates;
	recalibrate
} from "../software/remap/main"

workflow MAP
{
	take:
		genome_build
		ch_fasta
		ch_dir
		ch_common_vcf
		ch_data

	main:
		ch_data_branched = ch_data.map { it ->
			if ( it["type"] == "Normal" )
			{
				return (tuple (it, "Normal", it["normalR1"], it["normalR2"]))
			}
			if ( it["type"] == "Tumor" )
			{
				return (tuple (it, "Tumor", it["tumorR1"], it["tumorR2"]) )
			}
		}
		.dump (tag: 'remap input')
		.branch {
			paired: it[2].toString ().endsWith (".fastq.gz") && it[3].toString ().endsWith (".fastq.gz")
			other: true
		}

		ch_data_branched.other.view { "[MoCaSeq] error: Failed to find matching MAP workflow path for input:\n${it}" }

		fastqc_paired_extracted (genome_build, ch_data_branched.paired)
		trim_paired (fastqc_paired_extracted.out.result)

		trim_paired.out.result.set { ch_trim }

		fastqc_paired_trimmed (genome_build, ch_trim)
		bwa_mem_paired (ch_fasta, ch_trim)
		mark_duplicates (genome_build, bwa_mem_paired.out.result)
		recalibrate (genome_build, ch_fasta, ch_common_vcf, mark_duplicates.out.result)
		sample = recalibrate.out.result.map {
			it[0][it[0][[it[1].toLowerCase (),"BAM"].join ("")]] = it[2]
			it[0][it[0][[it[1].toLowerCase (),"BAI"].join ("")]] = it[3]
			it[0]
		}.reduce ( [:] ) { accumulator, item ->
			// Group by sample group
			if ( accumulator.containsKey (item["sampleGroup"]) )
			{
				accumulator[item["sampleGroup"]].add (item)
			}
			else
			{
				accumulator[item["sampleGroup"]] = [item]
			}
			accumulator
		}
		.flatMap ().flatMap { it ->
			// Within the sample group annotate Tumor samples with the new Normal bam
					def samples_by_type = it.value.inject ([:]) { accumulator, item ->
				if ( accumulator.containsKey (item["type"]) )
				{
					accumulator[item["type"]].add (item)
				}
				else
				{
					accumulator[item["type"]] = [item]
				}
				accumulator
			}
			println ("samples_by_type from map")
			println (samples_by_type)

			def result = []
			if ( samples_by_type.containsKey ("Normal") && samples_by_type.containsKey ("Tumor") )
			{
				result.addAll (samples_by_type["Normal"].collect { jt -> jt.putAll (["sampleGroupSize":it.value.size ()]);jt})
				result.addAll (samples_by_type["Tumor"].collect { jt -> jt.putAll (samples_by_type["Normal"][0].subMap (["normalBAM", "normalBAI"]));jt.put ("sampleGroupSize",it.value.size ());jt } )
			}
			else
			{
				samples_by_type.values ().each { jt -> result.addAll ( jt.collect { kt -> kt.put ("sampleGroupSize",it.value.size ());kt } ) }
			}
			result

		}.dump (tag: 'MAP output')

	emit:
		result = sample
}

workflow REMAP
{
	take:
		genome_build
		ch_fasta
		ch_dir
		ch_common_vcf
		ch_data

	main:
		ch_data_branched = ch_data.map { it ->
			if ( it["type"] == "Normal" )
			{
				return (tuple (it, "Normal", it["normalR1"], it["normalR2"]))
			}
			if ( it["type"] == "Tumor" )
			{
				return (tuple (it, "Tumor", it["tumorR1"], it["tumorR2"]) )
			}
		}
		.dump (tag: 'remap input')
		.branch {
			paired: it[2].toString ().endsWith (".bam") && it[2] == it[3]
			other: true
		}

		ch_data_branched.other.view { "[MoCaSeq] error: Failed to find matching REMAP workflow path for input:\n${it}" }

		sam_to_fastq_paired (ch_data_branched.paired.map { tuple (it[0], it[1], it[2] ) })
		fastqc_paired_extracted (genome_build, sam_to_fastq_paired.out.result)
		bwa_mem_paired (ch_fasta, fastqc_paired_extracted.out.result)
		mark_duplicates (genome_build, bwa_mem_paired.out.result)
		recalibrate (genome_build, ch_fasta, ch_common_vcf, mark_duplicates.out.result)
		sample = recalibrate.out.result.map { [it[0]["sampleName"], it] }
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

