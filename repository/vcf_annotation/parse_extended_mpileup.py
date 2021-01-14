#!/usr/local/bin/python
# coding: utf-8
"""
***********************************************
- PROGRAM: parse_extended_mpileup.py
- CONTACT: Gaurav Jain(gaurav.jain@tum.de)
***********************************************
"""
print(__doc__)

def main():

    # Get the indelbed dict
    indeldict = get_indelbed_dict(input_indelbed_file)

    # Get the pileupType
    ptype = 'Tumor'
    if pileupType =='n':
        ptype = 'Normal'

    # Import the extended pileup file
    print("\n- Importing the extended pileup file in chunks of 11 rows at a time...")
    try:
        extendedpileupDF = pd.read_csv(input_expileup_file, sep='\t', header=None, quoting=3)
        extendedpileupDF[0] = extendedpileupDF[0].astype(str)

        # head mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_mit_annotation_indels.extendedpileup | cut - f1-7
        # Final column is the ReadName
        # ┌───────┬────────────┬───────────────┬───────────────┬──────────────────────┬──────────────────────┬──────────────────────┐
        # │ Chrom │ GenomicPos │ Ref           │ TotalCov      │      ReadBases       │     BaseQ            │          MapQ        │
        # ├───────┼────────────┼───────────────┼───────────────┼──────────────────────┼──────────────────────┼──────────────────────┤
        # │     1 │  7411767   │ A             │            20 │ ,...,.......,....... │ ???>=>>>>>>;?>>>>=>< │ ]]]]]]]]]]]]]UUUUUUL │
        # │     1 │  7649718   │ c             │            15 │ ,.,.$A,,,,,,,a,,     │ BBBU@B@CCBBB?!B      │ ]SR@+J?E];]L4@V      │
        # └───────┴────────────┴───────────────┴───────────────┴──────────────────────┴──────────────────────┴──────────────────────┘
        # ...
        # Column 5: The bases string
        #     . (dot) means a base that matched the reference on the forward strand
        #     , (comma) means a base that matched the reference on the reverse strand
        #     </> (less-/greater-than sign) denotes a reference skip. This occurs, for example, if a base in the reference genome is intronic and a read maps to
        #         two flanking exons. If quality scores are given in a sixth column, they refer to the quality of the read and not the specific base.
        #     AGTCN denotes a base that did not match the reference on the forward strand
        #     agtcn denotes a base that did not match the reference on the reverse strand
        #     A sequence matching the regular expression \+[0-9]+[ACGTNacgtn]+ denotes an insertion of one or more bases starting from the next position
        #     A sequence matching the regular expression -[0-9]+[ACGTNacgtn]+ denotes a deletion of one or more bases starting from the next position
        #     ^ (caret) marks the start of a read segment and the ASCII of the character following `^' minus 33 gives the mapping quality
        #     $ (dollar) marks the end of a read segment
        #     * (asterisk) is a placeholder for a deleted base in a multiple basepair deletion that was mentioned in a previous line by the -[0-9]+[ACGTNacgtn]+ notation

        # List of Indels (tuples)
        indelList = []

        # Each chunk is in df format
        # Make a list of tuples with your data and then create a DataFrame with it:
        print("\n- Make a list of tuples with your data and then create a DataFrame with it:")
        print("\t- Get the chromsome and pos from the chunk of 11 rows")
        print("\t- Get the list of count of different INDEL types and list of different INDEL types")
        print("\t- Append the columns as the list elements into the final list")
        for k in indeldict.keys():
            chrom      = k.split(':')[0]
            pos        = k.split(':')[1]
            chunkStart = np.int(pos) - 5
            chunkEnd   = np.int(pos) + 5
            # print(chrom, pos, chunkStart, chunkEnd)
            chunk = extendedpileupDF.loc[(extendedpileupDF[0]=='{0}'.format(chrom)) & (extendedpileupDF[1] >= chunkStart) & (extendedpileupDF[1] <= chunkEnd)]

            # variantKey = 1:36706729
            variantKey = k

            # indelReadBases_Normal = ................|................|................|................^I.^].|................T.|........-2AC..........|.......*..........|.......*..........^?.|...................^].|.$...................|...................
            indelReadBases = "|".join(chunk[4])

            # indelPos_Normal = 36706724|36706725|36706726|36706727|36706728|36706729|36706730|36706731|36706732|36706733|36706734
            indelPos       = "|".join(chunk[1].map(str))

            # Get the list of count of different INDEL types and list of different INDEL types
            # itcln      = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
            # itln       = ['', '', '', '', '', '2AC', '', '', '', '', '']
            # ill        = [25, 28, 29, 36]
            # cleanBases = ',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,g,,,g'
            itcln,itln,_ill,_cleanBases = map(list, zip(*list(map(get_clean_readbase_indel_counts, list(chunk[4])))))
            indelTypeCountList        = "|".join(map(str,itcln))
            indelTypeList             = "|".join(itln)
            indelTypeCount            = np.nansum(itcln)

            # Append the columns as the list elements into the final list
            indelList.append((chrom, pos, variantKey, indelTypeCount, indelTypeCountList, indelTypeList, indelPos, indelReadBases))

        # Get the list of tuples into the dataframe
        print("\n- Get the list of tuples into the dataframe")
        parsedDF = pd.DataFrame(indelList, columns=('Chrom', 'GenomicPos', 'SearchID', '{0}ExtendedIndelTypeCount'.format(ptype), '{0}ExtendedIndelTypeCountList'.format(ptype), '{0}ExtendedIndelSubsType'.format(ptype), '{0}ExtendedIndelSubsPos'.format(ptype), '{0}ExtendedIndelSubsReadBases'.format(ptype)))
    except pd.errors.EmptyDataError:
        print("NOTE: {0} was empty. No usable indels .....Skipping....\n".format(input_expileup_file))
        # Create an empty output DF and add the relevant columns to it
        parsedDF = pd.DataFrame(columns=('Chrom', 'GenomicPos', 'SearchID', '{0}ExtendedIndelTypeCount'.format(ptype), '{0}ExtendedIndelTypeCountList'.format(ptype), '{0}ExtendedIndelSubsType'.format(ptype), '{0}ExtendedIndelSubsPos'.format(ptype), '{0}ExtendedIndelSubsReadBases'.format(ptype)))

    print("\n- Final shape of the dataframe: {}".format(parsedDF.shape))

    # Get output file name
    output_file = "{0}_parsed_extendedpileup.txt".format(get_file_info(input_expileup_file)[3])

    # Save the dataframe to output parsed file
    print("\n- Save the dataframe to output parsed file")
    parsedDF.to_csv(output_file, sep="\t", index=False, header=True,  float_format='%.f')

    print("\n- Your parsed output file is: {0}".format(output_file))

################ USER DEFINED FUNCTIONS ###################
def get_clean_readbase_indel_counts(bases):
    """
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
    indelList = list(chain.from_iterable(indelList))[1:]

    # This site has a 10 bp long Indel, but at the TumorExactIndelSubsType it is wrongly called as two Indels: -1|0tgtgtgtgta.
    # Correct would be: -10tgtgtgtgta
    # Location: 9_42008058_GTGTGTGTGTA_G
    # Merge the indexes correctly
    # input : ['+1a','-1', '0tgtgtgtgta','+2t', '-3','0ss']
    # output: ['+1a', '-10tgtgtgtgta', '+2t', '-30ss']
    idx    = []
    idxLoc = 0
    for i in indelList:
        if not i.islower():
            idx.append(idxLoc)
        idxLoc+=1
    xl = indelList.copy()

    # Merge the elements
    for j in idx:
        xl[j:j+1] = [''.join(xl[j:j+2])]

    # Remove the non-merged entries
    for i in sorted(idx, reverse=True):
        del xl[i+1]
    indelList = get_unique_list(xl.copy())

    # Calculate the count of indels only
    indelTypeCount = int(len(indelList))

    # Get substitutions
    # 16_30547590_GCACA_G
    # bases=',$,,,,,,,,g,,,,.,,,,-2ca,-2ca.,,,,t'
    # subsList = ','.join(map(str, re.findall('[acgt]', callString))) # ['g', 't']
    subsList = get_unique_list(re.findall('[acgt]', callString)) # ['g', 't']

    # Add the substituions along with the types of indels
    indelList.extend(subsList)

    # Get the length of the list to obtain number of unique indels
    # return(len(indelList), callString)
    return(indelTypeCount, "|".join(indelList), indelLocList, callString)

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

def get_indelbed_dict(input_indelbed_file):
    ''' Get the annotation in a dictionary '''

    # Get the information in a dictionary
    d = defaultdict(int)
    with open(input_indelbed_file, 'r') as fg:
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

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python scripts/parse_extended_mpileup.py -ip=output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2.extendedpileup -ib=output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_indels.bed
        -------------------------------------------------
        CONTACT:
            Gaurav Jain
            gaurav.jain@tum.de
        -------------------------------------------------
        '''))

    # Add arguments
    parser.add_argument("-ip", metavar='--ipfile', help="*Input extended pileup file", dest="input_expileup_file", type=str, required=True)
    parser.add_argument("-ib", metavar='--ibfile', help="*Input indelbed file       ", dest="input_indelbed_file", type=str, required=True)
    parser.add_argument('-pt',         '--prtype', help=" Enter pileup type file:\n\t -pt='t' for Tumor\n\t -pt='n' for Normal", dest="pileupType", default='t', action='store', choices=['t', 'n'])

    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().input_expileup_file:
        logdir="{0}/logs".format(get_file_info(get_file_info(parser.parse_args().input_expileup_file)[0])[0])
        create_dir(logdir)
        logfile = "{0}/{1}_parse_extended_mpileup.log".format(logdir, get_file_info(parser.parse_args().input_expileup_file)[1])
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
    import dask.dataframe as dd
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
    input_expileup_file = args.input_expileup_file
    input_indelbed_file = args.input_indelbed_file
    pileupType          = args.pileupType

    # Call the main function
    main()
