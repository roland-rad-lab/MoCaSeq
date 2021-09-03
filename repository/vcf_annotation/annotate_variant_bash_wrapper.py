#!/usr/local/bin/python
# coding: utf-8
"""
***********************************************
- PROGRAM: annotate_variant_bash_wrapper.py
***********************************************
"""

def main():

    # Set folder with all scripts
    scriptsRepo="/opt/MoCaSeq/repository/vcf_annotation/"

    # Get output file name
    scriptDir  = "{2}/scripts/wrapper/{1}".format(os.getcwd(), projName, projDir); os.system("mkdir -p {0}".format(scriptDir))
    scriptFile = "{0}/{1}_Mutect2_annotation_wrapper.sh".format(scriptDir, sampleName)
    ofile      = open(scriptFile, 'w')

    # Create a meta log file
    metaLogFile = "{0}/logs/{1}_annotation.log".format(projDir, sampleName); create_dir("{0}/logs".format(projDir))

    # Log everything in the meta log file
    # https://serverfault.com/questions/103501/how-can-i-fully-log-all-bash-scripts-actions
    # Explanation:
    # 1.  exec 3>&1 4>&2
    #     Saves file descriptors so they can be restored to whatever they were before redirection
    #     or used themselves to output to whatever they were before the following redirect.
    # 2.  trap 'exec 2>&4 1>&3' 0 1 2 3
    #     Restore file descriptors for particular signals. Not generally necessary since they should
    #     be restored when the sub-shell exits.
    # 3.  exec 1>log.out 2>&1
    #     Redirect stdout to file log.out then redirect stderr to stdout. Note that the order is
    #     important when you want them going to the same file. stdout must be redirected before stderr is redirected to stdout.

    ofile.write("#!/bin/bash\n\nexec 3>&1 4>&2\ntrap 'exec 2>&4 1>&3' 0 1 2 3\nexec 1> {} 2>&1\n".format(metaLogFile))
    # Everything below will go to the log file

    ofile.write("\n########################################\n# Log file: logs/{0}_annotation.log\npython -VV\n########################################\n".format(sampleName))
    # 1.1) Filter the input variant file for TumorVAFCustom == 0
    ofile.write("\n#- 1.1) Filter the input variant file for TumorVAFCustom == 0\n")
    ofile.write("echo '#- 1.1) Filter the input variant file for TumorVAFCustom == 0'\n")
    inputFilVariantFile     = "{0}_vafFiltered.{1}".format(interimFileName, get_file_info(inputVariantFile)[2])
    print("- inputFilVariantFile     = {0}".format(inputFilVariantFile))
    ofile.write("time python3.7 {2}/filter_variant_input.py -if={0} -of={1}\n".format(inputVariantFile, interimFileName, scriptsRepo))

    # Get the filteredInterimFileName
    filteredInterimFileName = "{0}_vafFiltered".format(interimFileName)
    print("- filteredInterimFileName = {0}".format(filteredInterimFileName))

    # 1.2) Convert merged vcf file to bed file
    ofile.write("\n#- 1.2) Convert merged vcf file to bed file\n")
    ofile.write("echo '#- 1.2) Convert merged vcf file to bed file'\n")
    snvBedFile     = "{}.bed".format(filteredInterimFileName)
    snvUcscBedFile = "{}.ucscbed".format(filteredInterimFileName)
    ofile.write("time python3.7 {2}/export_variant_input_to_bed.py -if={0} -of={1}\n".format(inputFilVariantFile, filteredInterimFileName, scriptsRepo))
    ofile.write("egrep -v 'Un|random|GL|JH|M|KQ|KZ|KI|EBV|Phix' {0} > {1} && mv {1} {0}\n".format(snvBedFile,"{0}.tmp".format(snvBedFile)))
    ofile.write("egrep -v 'Un|random|GL|JH|M|KQ|KZ|KI|EBV|Phix' {0} > {1} && mv {1} {0}\n".format(snvUcscBedFile,"{0}.tmp".format(snvUcscBedFile)))
    ofile.write("rm -rf {0} {1}\n".format("{0}.tmp".format(snvBedFile), "{0}.tmp".format(snvUcscBedFile)))

    # 1.3) Get the the subset bam for variant reads only
    ofile.write("\n#\t- 1.3) Get the the subset bam for variant reads only\n")
    ofile.write("echo '#\t- 1.3) Get the the subset bam for variant reads only'\n")
    snvSubBamFile = "{}_subset_variant.bam".format(filteredInterimFileName)
    ofile.write("time intersectBed -abam {0} -b {1} -wa > {2}\n".format(tumorBamFile, snvBedFile, snvSubBamFile))

    # 1.4) Make genome information (chromosome sizes) from the BAM file to get a genome.info file
    ofile.write("\n#- 1.4) Make genome information (chromosome sizes) from the BAM file to get a genome.info file\n")
    ofile.write("echo '#- 1.4) Make genome information (chromosome sizes) from the BAM file to get a genome.info file'\n")
    genomeInfoFile     = "{}_genome.info".format(filteredInterimFileName)
    genomeInfoUcscFile = "{}_genome.ucscinfo".format(filteredInterimFileName)
    ofile.write("samtools view -H {0}| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){{print $1,\"\\t\",$2,\"\\n\"}}' | egrep -v 'Un|random|GL|JH|M|KQ|KZ|KI|EBV|Phix' | sort -k1,1V -k2,2g -k3,3g > {1}\n".format(tumorBamFile, genomeInfoFile))
    ofile.write("samtools view -H {0}| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){{print \"chr\",$1,\"\\t\",$2,\"\\n\"}}' | egrep -v 'Un|random|GL|JH|M|KQ|KZ|KI|EBV|Phix' | sort -k1,1V -k2,2g -k3,3g > {1}\n".format(tumorBamFile, genomeInfoUcscFile))

    # 2) Get the base reads quality and mapping quality for the variants
    ofile.write("\n#- 2) Get the base reads quality and mapping quality for the variants\n")
    ofile.write("echo '#- 2) Get the base reads quality and mapping quality for the variants'\n")
    # 2.1) Get the pileup information for the variants
    ofile.write("\n#\t- 2.1) Get the pileup information for the variants\n")
    ofile.write("echo '#\t- 2.1) Get the pileup information for the variants'\n")
    snvMpileupFile  =  "{}.pileup".format(filteredInterimFileName)
    ofile.write("time samtools mpileup -A -B -Q 0 -s -x --output-QNAME -l {0} -f {1} {2} -o {3}".format(snvBedFile, genomeFastaFile, tumorBamFile,snvMpileupFile))
    # 2.2) Parse the output of exact tumor pileup
    ofile.write("\n\n#\t- 2.2) Parse the output of exact tumor pileup\n")
    ofile.write("echo '#\t- 2.2) Parse the output of exact tumor pileup'\n")
    snvPileupFile  = "{}_parsed_tumor_exact_pileup.txt".format(filteredInterimFileName)
    ofile.write("time python3.7 {2}/parse_mpileup.py -ip={0} -ib={1}\n".format(snvMpileupFile, snvBedFile, scriptsRepo))

    # 3) Segmental duplication
    ofile.write("\n#- 3) Get the segmental duplication information\n")
    ofile.write("echo '#- 3) Get the segmental duplication information'\n")
    segmentalDupBedFile = "{}_segmentalDups.bed".format(filteredInterimFileName)
    tempSegDupfile      = get_temp_file()
    ofile.write("egrep -v 'Un|random|GL|JH|M|KQ|KZ|KI|EBV|Phix' {0} > {1}\n".format(genomeSegDupFile, tempSegDupfile))
    ofile.write("time intersectBed -a {0} -b {1} -c -g {3}| sed 's/chr//' > {2}\n".format(snvUcscBedFile, tempSegDupfile, segmentalDupBedFile, genomeInfoUcscFile))
    ofile.write("\n#\t- Adding the header to the file\n")
    ofile.write("sed -i '1iChrom\\tStart\\tGenomicPos\\tVariantIDKey\\tSegmentalDupsCount' {}\n".format(segmentalDupBedFile))

    # 4) Repeat masker annotation
    ofile.write("\n#- 4) Get the repeat masker annotation\n")
    ofile.write("echo '#- 4) Get the repeat masker annotation'\n")
    repeatMaskerBedFile = "{}_repeatmasker.txt".format(filteredInterimFileName)
    ofile.write("time intersectBed -a {0} -b {1} -wao -g {3} > {2}\n".format(snvBedFile, genomeRepMaskFile, repeatMaskerBedFile, genomeInfoFile))
    ofile.write("\n#\t- Adding the header to the file\n")
    ofile.write("sed -i '1iChrom\\tStart\\tGenomicPos\\tVariantIDKey\\tRepChrom\\tRepGenoStart\\tRepGenoEnd\\tRepNameClassFamily\\tSwScore\\tSwStrand\\tRepMaskFlag' {}\n".format(repeatMaskerBedFile))

    # 5) GC percent annotation
    ofile.write("\n#- 5) GC percent annotation\n")
    ofile.write("echo '#- 5) GC percent annotation'\n")
    # 5.1) For snvBinGCpc
    ofile.write("\n#\t- 5.1) Get the GC% for the 5bp bin the SNV overlaps\n")
    ofile.write("echo '#\t- 5.1) Get the GC% for the 5bp bin the SNV overlaps'\n")
    snvBinGCpcFile = "{}_5bp_GC_percentage.txt".format(filteredInterimFileName)
    ofile.write("time intersectBed -a {0} -b {1} -wb -g {3} > {2}\n".format(genome5bpGCpcFile, snvUcscBedFile, snvBinGCpcFile, genomeInfoUcscFile))
    ofile.write("sed -i 's/chr//g' {}".format(snvBinGCpcFile))
    ofile.write("\n#\t- Adding the header to the file\n")
    ofile.write("sed -i '1iGCpc5bpChrom\\tGCpc5bpStart\\tGCpc5bpEnd\\tGCVariantBasePercentage\\tChrom\\tStart\\tGenomicPos\\tVariantIDKey' {}\n".format(snvBinGCpcFile))

    # 5.2) For snvReadMeanGCpc
    ofile.write("\n#\t- 5.2) Get the mean GC% for all the reads the SNV overlaps\n")
    ofile.write("#\t-      Map the reads to the snv locations, concatenate all the reads and calculate the mean gc% of all the mapped reads\n")
    ofile.write("echo '#\t- 5.2) Get the mean GC% for all the reads the SNV overlaps'\n")
    ofile.write("echo '#\t-      Map the reads to the snv locations, concatenate all the reads and calculate the mean gc% of all the mapped reads'\n")

    #snvReadMeanGCpcFile = "{}_snvReads_Mean_GC_percentage.txt".format(filteredInterimFileName)
    #ofile.write('time bedtools map -prec 2 -a {0} -b {1} -c 10,10 -o count,concat -g {3} | awk -v OFS="\\t" \'{{n=length($6); gc=gsub("[gcGC]", "", $6); print $1,$2,$3,$4,gc/n}}\'  > {2}\n\n'.format(snvBedFile, snvSubBamFile, snvReadMeanGCpcFile, genomeInfoFile))

    snvReadMeanGCpcFileRAW = "{}_snvReads_Mean_GC_percentage_raw.txt".format(filteredInterimFileName)
    ofile.write('time bedtools map -prec 2 -a {0} -b {1} -c 10,10 -o count,concat -g {3} | awk -v OFS="\\t" \'{{print $1,$2,$3,$4,$6}}\' > {2}\n\n'.format(snvBedFile, snvSubBamFile, snvReadMeanGCpcFileRAW, genomeInfoFile))

    snvReadMeanGCpcFileSTRING = "{}_snvReads_Mean_GC_percentage_string.txt".format(filteredInterimFileName)
    ofile.write('awk -v OFS="\\t" \'{{print $5}}\' {0} > {1}\n\n'.format(snvReadMeanGCpcFileRAW, snvReadMeanGCpcFileSTRING))

    snvReadMeanGCpcFilePERCENT = "{}_snvReads_Mean_GC_percentage_percent.txt".format(filteredInterimFileName)
    ofile.write('while read p; do len=$(echo $p | sed \'s/N//g\' | tr -d \'\\n\' | wc -c); cnt=$(echo $p | grep -oh \'C\|G\|g\|c\' | tr -d \'\\n\' | wc -c); gc=$(awk "BEGIN {{printf \\"%.6f\\",${{cnt}}/${{len}}}}"); echo -e $gc; done<{0} > {1}\n\n'.format(snvReadMeanGCpcFileSTRING, snvReadMeanGCpcFilePERCENT))

    snvReadMeanGCpcFile = "{}_snvReads_Mean_GC_percentage.txt".format(filteredInterimFileName)
    ofile.write('paste {0} {1} | awk \'{{print $1,$2,$3,$4,$6}}\' > {2}\n\n'.format(snvReadMeanGCpcFileRAW, snvReadMeanGCpcFilePERCENT, snvReadMeanGCpcFile))

    ofile.write("\n\n#\t- Adding the header to the file\n")
    ofile.write("sed -i '1iChrom\\tStart\\tGenomicPos\\tVariantIDKey\\tVariantReadsMeanGCPc' {}\n".format(snvReadMeanGCpcFile))

    # 6) Get the 10 bp sequence around the SNV from the genome
    ofile.write("\n#- 6) Get the 10 bp flanking sequence around the SNV from the genome\n")
    ofile.write("echo '#- 6) Get the 10 bp flanking sequence around the SNV from the genome'\n")
    tenBpAroundSNVFile  = "{}_10bp_around_SNV.txt".format(filteredInterimFileName)
    temp10bpbinBaseFile = get_temp_file()
    ofile.write("time awk -v OFS=\"\\t\" '{{print $1,$2-10,$3+10,$3}}' {0} > {1}\n".format(snvBedFile, temp10bpbinBaseFile))
    ofile.write("time bedtools getfasta -fi {0}  -bed {1} -name -bedOut > {2}\n".format(genomeFastaFile, temp10bpbinBaseFile, tenBpAroundSNVFile))
    ofile.write("\n#\t- Adding the header to the file\n")
    ofile.write("sed -i '1iChrom\\tStart\\tEnd\\tGenomicPos\\tVariant10bpFlankingSequence' {}\n".format(tenBpAroundSNVFile))

    # 7) Get the distance to the closest probe from the SNV
    # WESProbeAnnotation = wesProbe_probeID|TargetID|Sequence|Replication
    ofile.write("\n#- 7) Get the distance to the closest probe from the SNV\n")
    ofile.write("echo '#- 7) Get the distance to the closest probe from the SNV'\n")
    closestProbeFile = "{}_closest_wes_probe.txt".format(filteredInterimFileName)
    tempclosestprobe = get_temp_file()
    tmpsortingfile   = get_temp_file()
    ofile.write("grep -v ^M {0} > {1}\n".format(wesProbesFile, tempclosestprobe))
    ofile.write("(head -n 2 {0} && tail -n +3 {0} | sort -k1,1V -k2,2g -k3,3g) > {1} && mv {1} {0}\n".format(tempclosestprobe, tmpsortingfile))
    ofile.write("time bedtools closest -a {0} -b {1} -d -t first -g {3}| uniq > {2}\n".format(snvBedFile, tempclosestprobe, closestProbeFile, genomeInfoFile))
    ofile.write("\n#\t- Adding the header to the file\n")
    ofile.write("sed -i '1iChrom\\tStart\\tGenomicPos\\tVariantIDKey\\tWESProbeChrom\\tWESProbeStart\\tWESProbeEnd\\tWESProbeAnnotation\\tWESProbeScore\\tWESProbeStrand\\tWESProbeDistanceToVariant' {}\n".format(closestProbeFile))

    # 8.1) Make genome information (chromosome sizes) from the Normal BAM file to get a genome.info file
    ofile.write("\n#- 8.1) Make genome information (chromosome sizes) from the Normal BAM file to get a genome.info file'\n")
    ofile.write("echo '#- 8.1) Make genome information (chromosome sizes) from the Normal BAM file to get a genome.info file'\n")
    normGenomeInfoFile     = "{}_genome_normal.info".format(filteredInterimFileName)
    normGenomeInfoUcscFile = "{}_genome_normal.ucscinfo".format(filteredInterimFileName)
    ofile.write("samtools view -H {0}| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){{print $1,\"\\t\",$2,\"\\n\"}}' | egrep -v 'Un|random|GL|JH|M|KQ|KZ|KI|EBV|Phix' > {1}\n".format(normalBamFile, normGenomeInfoFile))
    ofile.write("samtools view -H {0}| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){{print \"chr\",$1,\"\\t\",$2,\"\\n\"}}' | egrep -v 'Un|random|GL|JH|M|KQ|KZ|KI|EBV|Phix' > {1}\n".format(normalBamFile, normGenomeInfoUcscFile))

    # 8.2) Get usable indels list
    ofile.write("\n#- 8.2) Get usable tumor indels information\n")
    ofile.write("echo '#- 8.2) Get usable tumor indels information'\n")
    tumorIndelFile     = "{}_indels.txt".format(filteredInterimFileName)
    usableIndelBed     = "{}_indels.bed".format(filteredInterimFileName)
    usableIndelUcscBed = "{}_indels.ucscbed".format(filteredInterimFileName)
    ofile.write("time samtools index {0}\n".format(snvSubBamFile))
    ofile.write("time python3.7 {4}/get_usable_indels_for_annoatation.py -if={0} -pf={1} -bm={2} -od={3}\n".format(inputFilVariantFile, snvPileupFile, snvSubBamFile, interimSampleDir, scriptsRepo))

    # 8.3) Increase the size of each feature in 5 bases on both sides
    ofile.write("\n#- 8.3) Extend each feature by 5 bases on both sides\n")
    ofile.write("echo '#- 8.3) Extend each feature by 5 bases on both sides'\n")
    usableIndelExtendedBed     = "{}_indels.extendedbed".format(filteredInterimFileName)
    usableIndelExtendedUcscBed = "{}_indels.extendeducscbed".format(filteredInterimFileName)
    ofile.write("time bedtools slop -i {0} -g {1} -b 5 > {2}\n".format(usableIndelBed, normGenomeInfoFile, usableIndelExtendedBed))
    ofile.write("time bedtools slop -i {0} -g {1} -b 5 > {2}\n".format(usableIndelUcscBed, normGenomeInfoUcscFile, usableIndelExtendedUcscBed))

    # 9) Get the extended base reads and mapping quality for the variants for indels in tumor and normal bam files
    ofile.write("\n#- 9) Get the extended base reads and mapping quality for the variants for indels in tumor and normal bam files'\n")
    ofile.write("echo '#- 9) Get the extended base reads and mapping quality for the variants for indels in tumor and normal bam files'\n")

    # 9.1) Get the extended pileup information for the variants in tumor bam file
    ofile.write("\n#- 9.1) Get the extended pileup information for the variants in tumor bam file\n")
    ofile.write("echo '#- 9.1) Get the extended pileup information for the variants in tumor bam file'\n")
    snvMextendedTumorPileupFile  =  "{}_tumor.extendedpileup".format(filteredInterimFileName)
    ofile.write("time samtools mpileup -A -B -Q 0 -s -x --output-QNAME -l {0} -f {1} {2} -o {3}".format(usableIndelExtendedBed, genomeFastaFile, tumorBamFile,snvMextendedTumorPileupFile))

    # 9.2) Parse the output of extended pileup file for tumor bams
    ofile.write("\n\n#- 9.2) Parse the output of extended pileup file for tumor bams\n")
    ofile.write("echo '#- 9.2) Parse the output of extended pileup file for tumor bams'\n")
    extParsedTumorPileupFile  = "{}_tumor_parsed_extendedpileup.txt".format(filteredInterimFileName)
    ofile.write("time python3.7 {2}/parse_extended_mpileup.py -ip={0} -ib={1} -pt=t \n".format(snvMextendedTumorPileupFile, usableIndelBed, scriptsRepo))

    # 9.3) Get the extended base reads and mapping quality for the variants for indels in normal bam files
    ofile.write("\n#- 9.3) Get the extended base reads and mapping quality for the variants for indels in normal bam files")
    ofile.write("echo '#- 9.3) Get the extended base reads and mapping quality for the variants for indels in normal bam files'\n")
    snvMextendedNormalPileupFile  =  "{}_normal.extendedpileup".format(filteredInterimFileName)
    ofile.write("time samtools mpileup -A -B -Q 0 -s -x --output-QNAME -l {0} -f {1} {2} -o {3}".format(usableIndelExtendedBed, genomeFastaFile, normalBamFile,snvMextendedNormalPileupFile))

    # 9.4) Parse the output of extended pileup file for normal bams
    ofile.write("\n\n#- 9.4) Parse the output of extended pileup file for normal bams")
    ofile.write("echo '#- 9.4) Parse the output of extended pileup file for normal bams'\n")
    extParsedNormalPileupFile  = "{}_normal_parsed_extendedpileup.txt".format(filteredInterimFileName)
    ofile.write("time python3.7 {2}/parse_extended_mpileup.py -ip={0} -ib={1} -pt=n \n".format(snvMextendedNormalPileupFile, usableIndelBed, scriptsRepo))

    # 10) Merge nnotations to the original file (inputVariantFile)
    ofile.write("\n#- 10) Merge annotations to the original file (inputVariantFile)\n")
    ofile.write("echo '#- 10) Merge annotations to the original file (inputVariantFile)'\n")
    ofile.write("time python3.7 {12}/merge_annotation_files.py -if={0} -pf={1} -sf={2} -rf={3} -gb={4} -ab={5} -gr={6} -cp={7} -ti={8} -tx={9} -nx={10} -sn={11}\n".format(inputVariantFile, snvPileupFile, segmentalDupBedFile, repeatMaskerBedFile, snvBinGCpcFile, tenBpAroundSNVFile, snvReadMeanGCpcFile, closestProbeFile, tumorIndelFile, extParsedTumorPileupFile, extParsedNormalPileupFile, sampleName, scriptsRepo))
    ofile.close()
    print("- Your script file is     = {0}".format(scriptFile))

################ USER DEFINED FUNCTIONS ###################

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        python3.7 annotate_variant_bash_wrapper.py -gf=/home/rad/users/gaurav/projects/pdacMetastasis/input/annotation/mm10/GRCm38.p6.fna -tb=/media/rad/HDD1/pdacMetastasis/input/mPDAC/cherrypicked1/DS08_5320_LivMet-1/results/bam/DS08_5320_LivMet-1.Tumor.bam -nb=/media/rad/HDD1/pdacMetastasis/input/mPDAC/cherrypicked1/DS08_5320_LivMet-1/results/bam/DS08_5320_LivMet-1.Normal.bam -sd=/home/rad/users/gaurav/projects/pdacMetastasis/input/annotation/mm10/mm10_ucsc_segmentalDups.bed -rp=/home/rad/users/gaurav/projects/pdacMetastasis/input/annotation/mm10/mm10_repeatMasker.bed -gc=/home/rad/users/gaurav/projects/pdacMetastasis/input/annotation/mm10/mm10_5bp_GC_percentage.bedgraph -pb=/home/rad/users/gaurav/projects/pdacMetastasis/input/annotation/mm10/wes_S0276129_Probes_mm10.bed -pn=mPDAC_cherrypicked1 -pd=/home/rad/users/gaurav/projects/pdacMetastasis/output/mPDAC_cherrypicked1
        -------------------------------------------------
        -------------------------------------------------
        '''))

    # Add arguments
    parser.add_argument("-pd", metavar='--prjdir', help="project output directory"    , dest="projDir"              , type=str, required=True)
    parser.add_argument("-sn", metavar='--sample', help="Name of the sample. ex: DS4" , dest="sampleName"           , type=str, required=True)
    parser.add_argument("-gf", metavar='--gffile', help="input genome fasta file"     , dest="genomeFastaFile"      , type=str, required=True)
    parser.add_argument("-tb", metavar='--tbfile', help="input tumor  bam file"       , dest="tumorBamFile"         , type=str, required=True)
    parser.add_argument("-nb", metavar='--nbfile', help="input normal bam file"       , dest="normalBamFile"        , type=str, required=True)
    parser.add_argument("-pn", metavar='--prname', help="project or experiment name"  , dest="projName"             , type=str, required=True)
    parser.add_argument("-sd", metavar='--sdfile', help="input segmental duplication file", dest="genomeSegDupFile" , type=str, required=True)
    parser.add_argument("-rp", metavar='--rpfile', help="input repeat masker file"    , dest="genomeRepMaskFile"    , type=str, required=True)
    parser.add_argument("-gc", metavar='--gcfile', help="input genomewide GC percent in 5-base windows bedgraph file", dest="genome5bpGCpcFile", type=str, required=True)
    parser.add_argument("-pb", metavar='--pbfile', help="Input probes bed file", dest="wesProbesFile", type=str, required=True)

    # http://thomas-cokelaer.info/blog/2014/03/python-argparse-issues-with-the-help-argument-typeerror-o-format-a-number-is-required-not-dict/

    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().projDir:
        logdir="{0}/interimFiles/logs".format(parser.parse_args().projDir)
        create_dir(logdir)
        logfile = "{0}/{1}_annotate_variant_wrapper.log".format(logdir, parser.parse_args().sampleName)
    else:
        logdir  = "{0}/logs".format(os.getcwd())
        create_dir(logdir)
        logfile = "{0}/{1}.log".format(logdir,get_file_info(sys.argv[0])[1])

    logf = open(logfile, 'w')
    sys.stdout = Log(logf, sys.stdout)

    # Parse command line with parse_args and store it in an object
    args = parser.parse_args()
    print_initial_arguments(parser)
    return args

if __name__=="__main__":
    from functools import partial, reduce
    print (__doc__)

    # Built in modules
    import argparse
    import os.path
    import sys

    # 3rd party modules
    import textwrap
    import re
    import numpy as np
    import scipy as sp
    import pandas as pd
    import matplotlib as mp
    #mp.use('Agg') # to use matplotlib without X11
    import matplotlib.pyplot as plt
    import subprocess
    import binascii as bi
    import scipy.stats as stats
    from collections import *
    from numpy import nanmean

    # for looping files in a dir
    import glob

    # user defined modules
    from gjainPyLib import *     # import all the functions from the Gaurav`s python scripts/library

    ### for color scale
    from  matplotlib import colors
    from itertools import cycle, islice # barplot colors

    ################ USER CONFIGURATION ###################
    np.set_printoptions(precision=6)
    #######################################################

    # Get input options
    args = check_options()

    # Store the variables
    genomeFastaFile   = args.genomeFastaFile
    genomeSegDupFile  = args.genomeSegDupFile
    genomeRepMaskFile = args.genomeRepMaskFile
    genome5bpGCpcFile = args.genome5bpGCpcFile
    wesProbesFile     = args.wesProbesFile
    tumorBamFile      = args.tumorBamFile
    normalBamFile     = args.normalBamFile
    projDir           = args.projDir
    projName          = args.projName
    sampleName        = args.sampleName

    # User variables
    interimSampleDir  = "{0}/interimFiles/{1}".format(projDir,sampleName); create_dir(interimSampleDir)
    interimFileName   = "{0}/{1}_Mutect2".format(interimSampleDir,sampleName)
    inputVariantFile  = "{0}/{1}_Mutect2.txt".format(projDir,sampleName)

    main()
