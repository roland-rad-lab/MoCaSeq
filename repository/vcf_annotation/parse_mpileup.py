#!/usr/local/bin/python
# coding: utf-8
"""
***********************************************
- PROGRAM: parse_mpileup.py
- CONTACT: Gaurav Jain(gaurav.jain@tum.de)
***********************************************
"""

def main():

    # Get the snvbed dict
    snvdict = get_snvbed_dict(input_snvbed_file)

    print("- Reading and converting mergedvcf file to bed file...")
    with open(input_pileup_file, 'r') as f:
        # http://samtools.sourceforge.net/pileup.shtml
        # head mergedvcf/txt/DS08_5320_LivMet-1-Mutect2.pileup | cut - f1-7
        # Final column is the ReadName
        # ┌───────┬──────────┬───────────────┬───────────────┬──────────────────────┬──────────────────────┬──────────────────────┐
        # │ Chrom │ Position │ ReferenceBase │ totalCoverage │      ReadBases       │     BaseQuality      │    MappingQuality    │
        # ├───────┼──────────┼───────────────┼───────────────┼──────────────────────┼──────────────────────┼──────────────────────┤
        # │     1 │  7411767 │ A             │            20 │ ,...,.......,....... │ ???>=>>>>>>;?>>>>=>< │ ]]]]]]]]]]]]]UUUUUUL │
        # │     1 │  7649718 │ c             │            15 │ ,.,.$A,,,,,,,a,,     │ BBBU@B@CCBBB?!B      │ ]SR@+J?E];]L4@V      │
        # └───────┴──────────┴───────────────┴───────────────┴──────────────────────┴──────────────────────┴──────────────────────┘
        # ...

        # Write header to the output line
        pileupHeader         = "Chrom\tGenomicPos\tRef\tTotalCov\tReadBases\tBaseQ\tMapQ\tReadName"
        BaseQualHeader       = "\t".join(["BaseQ" + l for l in ['Mean', 'Std', 'Min', '25pc', 'Median', '75pc', 'Max']])
        mappingQualityHeader = "\t".join(["MapQ" + l for l in ['Mean', 'Std', 'Min', '25pc', 'Median', '75pc', 'Max']])
        ofile.write("{0}\t{1}\t{2}\n".format(pileupHeader, BaseQualHeader, mappingQualityHeader))

        for line in useful_lines(f):
            # Split the line
            row  = line.split("\t")
            # Get relevant keys from pileup file
            pileupKey      = "{0}:{1}".format(row[0], row[1])  # 1:7649718
            _referenceBase = row[2].upper()                    # C
            _totalCoverage = row[3]                            # 15
            readBases      = row[4]                            # ,.,.$A,,,,,,,A,,
            baseQuality    = row[5]                            # BBBU@B@CCBBB?!B
            mappingQuality = row[6]                            # ]SR@+J?E];]L4@V
            itcln,_itln,ill,callString = get_clean_readbase_indel_counts(readBases)
            callString     = callString.upper()
            subsList       = get_unique_list(re.findall('[ACGT]', callString)) # ['g', 't']

            # Get the snvbed dictionary value for pileupKey
            snvbed = snvdict[pileupKey]                   # new: 1_3160218_T_C_4 ...old: DS08_5320_LivMet-1|1_7649718_CC_AT|2|DS08_5320_LivMet-3-Mutect2

            # Parse snvbed value to respective variables (m2  = mutect2)
            try:
                m2row            = snvbed.split('_')
                _m2ReferenceBase = m2row[2].upper()     # CC: refrence base ... store only first base (C) as second is not considered
                m2AltBase        = m2row[3].upper()     # AT: Alt      base ... store only first base (A) as second is not considered
                _m2AltReadsCount = m2row[4].upper()     # 4 : Alt reads (skip if zero)

                ## Get the index of the snv from the callString and then use that index to get the quality of those reads
                # Get the index of snv
                # NOTE: Before calculating the indexes, I need to remove $ and ^ sign which are the part of the CALF format
                #       Also at the read base column, a symbol `^' marks the start of a read segment which is a contiguous
                #       subsequence on the read separated by `N/S/H' CIGAR operations. The ASCII of the character following
                #       `^' minus 33 gives the mapping quality. A symbol `$' marks the end of a read segment. Start and end
                #       markers of a read are largely inspired by Phil Green's CALF (compact alignment format) format.
                #       These markers make it possible to reconstruct the read sequences from pileup.
                # http://www.phrap.org/phredphrap/calf.pdf

                # https://stackoverflow.com/questions/14215338/python-remove-multiple-character-in-list-of-string/14215379
                # Before: ,.,.$A,,,,,,,A,,
                # After : ,.,.A,,,,,,,A,,
                # Positions of A in ,.,.A,,,,,,,A,, are [4, 12]
                snvIdx = [pos for pos, c in enumerate(callString) if c == m2AltBase[0]]

                # debugPos = '7:119640568'
                # if pileupKey == debugPos:
                #     print(row)
                #     print("subsList      : {0}".format(subsList))
                # #     print("m2row         : {0}". format(m2row))
                #     print("snvIdx        : {0}". format(snvIdx))
                # #     print("m2AltBase     : {0}". format(m2AltBase))
                # #     print("readBases     : {0}". format(readBases))
                # #     print("callString    : {0}". format(callString))
                # #     print("mappingQuality: {0}". format(mappingQuality))
                # #     print("indelListCount: {0}". format(itcln))
                # #     print("indelList     : {0}". format(_itln))
                # #     print("indelLocation : {0}". format(ill))
                # #     print("callString    : {0}". format(callString))


                baseQscores          = "\t\t\t\t\t\t"
                mappingQualityscores = "\t\t\t\t\t\t"
                if len(subsList) > 0 and m2AltBase[0] in subsList:
                    # if pileupKey == debugPos:   print("1")
                    # if pileupKey == debugPos: print("here: {}".format(subsList))
                    # Get the base and mapping quality for the snv index
                    # https://www.geeksforgeeks.org/python-accessing-all-elements-at-given-list-of-indexes/
                    # https://stackoverflow.com/questions/25082410/apply-function-to-each-element-of-a-list?rq=1
                    # Decided to go for descripting of all the reads for that variant
                    # In [87]: pd.DataFrame([ord(i)-33 for i in list(itemgetter(*snvIdx)(mappingQuality))]).describe().T
                    #        count  mean      std       min    25%    50%   75%    max
                    #   0    2.0    14.5      6.363961  10.0   12.25  14.5  16.75  19.0
                    # Final string = '30.50\t0.71\t30\t30.25\t30.50\t30.75\t31'
                    baseQDF = pd.DataFrame([ord(i)-33 for i in list(itemgetter(*snvIdx)(baseQuality))]).describe().T.values.tolist()[0]
                    baseQscores = "{0:1.2f}\t{1:1.2f}\t{2:1.0f}\t{3:1.2f}\t{4:1.2f}\t{5:1.2f}\t{6:1.0f}".format(baseQDF[1], baseQDF[2], baseQDF[3], baseQDF[4], baseQDF[5], baseQDF[6], baseQDF[7])
                    mappingQualityDF  = pd.DataFrame([ord(i)-33 for i in list(itemgetter(*snvIdx)(mappingQuality))]).describe().T.values.tolist()[0]
                    # if pileupKey == debugPos:   print(mappingQuality)
                    mappingQualityscores  = "{0:1.2f}\t{1:1.2f}\t{2:1.0f}\t{3:1.2f}\t{4:1.2f}\t{5:1.2f}\t{6:1.0f}".format(mappingQualityDF[1], mappingQualityDF[2], mappingQualityDF[3], mappingQualityDF[4], mappingQualityDF[5], mappingQualityDF[6], mappingQualityDF[7])
                elif itcln > 0:
                    # if pileupKey == debugPos:   print("INDELs")
                    mappingQualityDF      = pd.DataFrame([ord(i)-33 for i in list(itemgetter(*ill)(mappingQuality))]).describe().T.values.tolist()[0]
                    mappingQualityscores  = "{0:1.2f}\t{1:1.2f}\t{2:1.0f}\t{3:1.2f}\t{4:1.2f}\t{5:1.2f}\t{6:1.0f}".format(mappingQualityDF[1], mappingQualityDF[2], mappingQualityDF[3], mappingQualityDF[4], mappingQualityDF[5], mappingQualityDF[6], mappingQualityDF[7])
                else:
                    # if pileupKey == debugPos:  print("{2}: Alteration {0} is not found in list of all alterations {1}".format(m2AltBase[0], subsList, pileupKey))
                    # print("{}".format(pileupKey))
                    print("{2}: Alteration {0} is not found in list of all alterations {1}".format(m2AltBase[0], subsList, pileupKey))
                    pass

                # Write to the output line
                ofile.write("{0}\t{1}\t{2}\n".format(line, baseQscores, mappingQualityscores))
            except:
                # print("MY INDEL ERROR:", line)
                pass

    ofile.close()
    print("\n- Your text  output file is: {0}".format(output_file))

################ USER DEFINED FUNCTIONS ###################
def get_snvbed_dict(input_snvbed_file):
    ''' Get the annotation in a dictionary '''

    # Get the information in a dictionary
    d = defaultdict(int)
    with open(input_snvbed_file, 'r') as fg:
        ######### OLD ##############
        # sort - V - k1, 1 - k3, 3 mergedvcf/txt/DS08_5320_LivMet-1-Mutect2.bed | head
        # ┌───┬──────────┬──────────┬─────────────────────────────────────────────────────────────────┐
        # │ 1 │  7411766 │  7411767 │ DS08_5320_LivMet-1|1_7411767_A_G|0|DS08_5320_LungMet-1-Mutect2  │
        # │ 1 │  7649717 │  7649718 │ DS08_5320_LivMet-1|1_7649718_CC_AT|2|DS08_5320_LivMet-3-Mutect2 │
        # │ 1 │ 10554584 │ 10554585 │ DS08_5320_LivMet-1|1_10554585_C_G|0|DS08_5320_LungMet-1-Mutect2 │
        # │ 1 │ 19228312 │ 19228313 │ DS08_5320_LivMet-1|1_19228313_AG_A|0|DS08_5320_LivMet-3-Mutect2 │
        # └───┴──────────┴──────────┴─────────────────────────────────────────────────────────────────┘
        # ....

        # NEW
        # ┌───┬─────────┬─────────┬──────────────────┐
        # │ 1 │ 3160217 │ 3160218 │ 1_3160218_T_C_4  │
        # │ 1 │ 3409856 │ 3409857 │ 1_3409857_A_G_22 │
        # │ 1 │ 5059015 │ 5059016 │ 1_5059016_T_A_11 │
        # └───┴─────────┴─────────┴──────────────────┘
        # ...


        # Loop through rest of the file
        for line in useful_lines(fg):
            row = re.split('\t', line)
            try:
                d["{0}:{1}".format(row[0], row[2])] = row[3]  # key = 1:7649718
            except:
                pass

    return d

def get_clean_readbase_indel_counts(bases):
    """
    Removes from the call string in,$,,,,,,,,,,,,,,,,,,,,,,,,,,-2tg,,,-2tg,-2tg,,,g,,,-2tgg mpileup (5th column) the ^ character and
    the char next to it.,$,,,,,,,,,,,,,,,,,,,,,,,,,,-2tg,,,-2tg,-2tg,,,g,,,-2tgg
    bases:,$,,,,,,,,,,,,,,,,,,,,,,,,,,-2tg,,,-2tg,-2tg,,,g,,,-2tgg
        String of read bases (5th col,$,,,,,,,,,,,,,,,,,,,,,,,,,,-2tg,,,-2tg,-2tg,,,g,,,-2tggumn of mpileup)
    Return:,$,,,,,,,,,,,,,,,,,,,,,,,,,,-2tg,,,-2tg,-2tg,,,g,,,-2tgg
        Same string as bases but with,$,,,,,,,,,,,,,,,,,,,,,,,,,,-2tg,,,-2tg,-2tg,,,g,,,-2tgg ^ and the following char removed as well as
        indels,$,,,,,,,,,,,,,,,,,,,,,,,,,,-2tg,,,-2tg,-2tg,,,g,,,-2tgg
    Example:,$,,,,,,,,,,,,,,,,,,,,,,,,,,-2tg,,,-2tg,-2tg,,,g,,,-2tgg
        bases= '^A....,,.,.,...,,,.,....^k.'
        get_clean_readbase(bases) >>> '....,,.,.,...,,,.,.....'
    Source: Dario Beraldi (https://github.com/dariober)
            https://raw.githubusercontent.com/dariober/bioinformatics-cafe/master/bam2methylation/bam2methylation.py

    Another good read: https://www.biostars.org/p/76829/

    Also:
    Count occurance of differnt types of indels in the base string
    For example:
    1) ,$,,,,,,,,,,,,,,,,,,,,,,,,,,-2tg,,,-2tg,-2tg,,,g,,,-2tgg                        has 1 types [-2tg]
    2) ,$,,,,,,,,,,,,,.,,.,,,,,,,,,,aa,,,,-4gaga,-2ga,,,,,,,-2ga,,,-4gaga,,,,-2gaa*,,, has 2 types ['4gaga', '2ga']

    """
    # NOTE:
    callString     = ''
    skip           = False
    getIndelNum    = False # Switch to start accumulating ints following +/- if indelFound flag is set
    indelNumList   = []    # List of ints following +/-. Converted to int() will give indel length
    nskip          = 0
    indelList      = []    # List of all indels (unique)

    # # Retun NAN if input is nan
    # if(pd.isna(bases)):
    #     print("\n- Skipping NANs ...")
    #     return (np.nan, np.nan)

    indelType    = ''  # + or -
    indelLoc     = 0
    indelLocList = []
    for i, x in enumerate(bases):
        # if bases == ',,,,,,+2tg,+2tg,,+2tg': print("\n{}) origLoc {}: {}".format(i, indelLoc, x))
        # print("\n{}) origLoc {}: {}".format(i, indelLoc, x))
        if nskip > 0:
            nskip -= 1
            indelList.append(x.lower())
            indelLoc -= 1
        elif x  == '^':
            skip= True
            indelLoc -= 1
        elif skip:
            skip= False
            indelLoc -= 1
        elif x == '$':
            indelLoc -= 1
            continue # Skip end-of-read marker
        elif x in ('+', '-'):
            indelLoc -= 1
            indelType = x
            indelLocList.append(indelLoc)
            getIndelNum = True
        elif getIndelNum:
            indelLoc -= 1
            if x.isdigit():
                indelNumList.append(x)
                indelList.append('|')
                indelList.append(indelType)
                indelList.append(x.lower())
                indelType=''
            else:
                indelList.append(x.lower())
                try:
                    nskip = int(''.join(indelNumList))-1
                except:
                    pass
                indelNumList= []
                getIndelNum= False
        else:
            # if bases == ',,,,,,+2tg,+2tg,,+2tg': print("{}) Location{}: {}".format(i, indelLoc, x))
            # print("{}) Location{}: {}".format(i, indelLoc, x))
            callString += x

        indelLoc += 1

    # Split the list into sublists separated by '|'
    # Input  = [',', '|', '2', 't', 'g', '|', '2', 't', 'g', '|', '2', 't', 'g', '|', '2', 't', 'g']
    # Output = [[','], ['2tg'], ['2tg'], ['2tg'], ['2tg']]
    # Or
    # Input  = [',', '|', '4', 'g', 'a', 'g', 'a', '|', '2', 'g', 'a', '|', '2', 'g', 'a', '|', '4', 'g', 'a', 'g', 'a', '|', '2', 'g', 'a']
    # Output = [[','], ['4gaga'], ['2ga'], ['2ga'], ['4gaga'], ['2ga']]
    indelList = [[''.join(list(g))] for k,g in groupby(breakby(indelList, '|', None), lambda x: x is not None) if k]

    # Flatten the list and remove the first element (',')
    # Input  = [[','], ['2tg'], ['2tg'], ['2tg'], ['2tg']]
    # Output = ['2tg']
    # Or
    # Input  = [[','], ['4gaga'], ['2ga'], ['2ga'], ['4gaga'], ['2ga']]
    # Output = ['4gaga', '2ga']
    indelList = get_unique_list(list(chain.from_iterable(indelList))[1:])

    # Get the length of the list to obtain number of unique indels
    # return(len(indelList), callString)
    return(int(len(indelList)), "|".join(indelList), indelLocList, callString)

def breakby(biglist, sep, delim=None):
    '''
    Split list into lists based on a character occurring inside of an element
    Source: https://stackoverflow.com/questions/45281189/split-list-into-lists-based-on-a-character-occurring-inside-of-an-element/45281271
    '''
    from itertools import groupby
    for item in biglist:
        p = item.split(sep)
        yield p[0]
        if len(p) > 1:
            yield delim
            yield p[1]

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python scripts/parse_mpileup.py -ip=output/mPDAC_cherrypicked1/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2.pileup -ib=output/mPDAC_cherrypicked1/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2.bed
        -------------------------------------------------
        CONTACT:
            Gaurav Jain
            gaurav.jain@tum.de
        -------------------------------------------------
        '''))

    # Add arguments
    parser.add_argument("-ip", metavar='--ipfile', help="input  pileup file", dest="input_pileup_file", type=str, required=True)
    parser.add_argument("-ib", metavar='--ibfile', help="input  snvbed file", dest="input_snvbed_file", type=str, required=True)

    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().input_pileup_file:
        logdir="{0}/logs".format(get_file_info(get_file_info(parser.parse_args().input_pileup_file)[0])[0])
        create_dir(logdir)
        logfile = "{0}/{1}_parse_tumor_exact_pileup.log".format(logdir, get_file_info(parser.parse_args().input_pileup_file)[1])
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
    print(__doc__)

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
    from operator import itemgetter
    from itertools import groupby, chain

    # for looping files in a dir
    import glob

    # user defined modules
    from gjainPyLib import *      # import all the functions from the Gaurav`s python library

    ### for color scale
    from  matplotlib import colors
    from itertools import cycle, islice # barplot colors

    ################ USER CONFIGURATION ###################
    np.set_printoptions(precision=6)
    #######################################################

    # Get input options
    args = check_options()

    # Store the variables
    input_pileup_file = args.input_pileup_file
    input_snvbed_file = args.input_snvbed_file

    # Get output file name
    output_file = "{0}_parsed_tumor_exact_pileup.txt".format(get_file_info(input_pileup_file)[3])
    ofile = open(output_file, 'w')

    main()
