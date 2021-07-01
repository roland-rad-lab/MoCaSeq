
include { jabba_matched } from "../software/jabba/main"

workflow JABBA
{
	take:
		ch_manta


	main:

		jabba_matched (ch_manta)
}


