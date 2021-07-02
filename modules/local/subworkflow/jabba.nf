
include { jabba_matched } from "../software/jabba/main"

workflow JABBA
{
	take:
		ch_manta
		ch_ratio

	main:

		ch_manta_and_ratio = ch_ratio.dumps (tag: "ch_ratio")

		jabba_matched (ch_manta)
}


