#!/usr/local/bin/python
# coding: utf-8
"""
***********************************************
- PROGRAM: export_variant_input_to_bed.py
***********************************************
"""

def main():

    # 1) Import original mergedvcf file
    # mergedvcfFile = 'output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2.txt'
    # ┌────────────────┬───────────────────────┬──────────────────────┬───────────────────────┬──────┐
    # │   Old_Header   │      New_Header       │         SNV1         │         SNV2          │ SNV… │
    # ├────────────────┼───────────────────────┼──────────────────────┼───────────────────────┼──────┤
    # │ fileID         │ FileID                │ DS08_5320_LivMet-1   │ DS08_5320_LivMet-1    │ …    │
    # │ ID             │ VariantID             │ 1_46083693_T_G       │ 1_66632104_C_T        │ …    │
    # │ SearchID       │ SearchID              │ 1:46083693           │ 1:66632104            │ …    │
    # │ CHROM          │ Chrom                 │ 1                    │ 1                     │ …    │
    # │ POS            │ GenomicPos            │ 46083693             │ 66632104              │ …    │
    # │ REF            │ Ref                   │ T                    │ C                     │ …    │
    # │ ALT            │ Alt                   │ G                    │ T                     │ …    │
    # │ GEN[Tumor].VF  │ TumorVAFMutect2       │ 0.129032             │ 0.2                   │ …    │
    # │ VAF            │ TumorVAFCustom        │ 0.129032258064516    │ 0.2                   │ …    │
    # │ totalCov       │ TotalCov              │ 62                   │ 20                    │ …    │
    # │ GEN[Tumor].RD  │ TumorRefCov           │ 54                   │ 16                    │ …    │
    # │ GEN[Tumor].AD  │ TumorAltCov           │ 8                    │ 4                     │ …    │
    # │ GEN[Normal].RD │ NormalRefCov          │ 97                   │ 16                    │ …    │
    # │ GEN[Normal].AD │ NormalAltCov          │ 0                    │ 0                     │ …    │
    # │ GENE           │ GeneName              │ Dnah7b               │ Unc80                 │ …    │
    # │ EFFECT         │ VariantLocationInGene │ intron_variant       │ intron_variant        │ …    │
    # │ FEATUREID      │ GeneID                │ ENSMUST00000069293.8 │ ENSMUST00000061620.15 │ …    │
    # │ HGVS_C         │ HGVSc                 │ c.274-12T>G          │ c.5812+573C>T         │ …    │
    # │ HGVS_P         │ HGVSp                 │                      │                       │ …    │
    # │ RESCUED        │ Rescued               │ no                   │ no                    │ …    │
    # └────────────────┴───────────────────────┴──────────────────────┴───────────────────────┴──────┘
    print("- Reading and converting mergedvcf file to bed file...")
    mergedvcfDF = pd.read_csv(inputVariantFile, sep='\t', quoting=3)

    # Create an empty output DF and add the relevant columns to it
    bedDF = pd.DataFrame()
    bedDF['Chrom'] = mergedvcfDF['Chrom']
    bedDF['start'] = mergedvcfDF['GenomicPos'] - 1
    bedDF['end']   = mergedvcfDF['GenomicPos']
    bedDF['chr_pos_ref_alt']  = mergedvcfDF.apply(lambda row: "{0}_{1:.0f}_{2}_{3}_{4:.0f}".format(row['Chrom'], row['GenomicPos'], row['Ref'], row['Alt'], row['TumorRefCov']), axis=1)

    # Save the dataframe to output bed file (Ensembl format without chr)
    print("\n- Saving the results to the output files:")
    bedDF.to_csv(output_file, sep="\t", index=False, header=False,  float_format='%.f')

    # Save the dataframe to ucscoutput bed file (UCSC bed format with chr)
    bedDF['Chrom'] = 'chr' + bedDF['Chrom'].astype(str)
    bedDF.to_csv(ucscoutput_file, sep="\t", index=False, header=False,  float_format='%.f')

    # Sort the bed files
    print("\n- Sort the bed files")
    os.system("sort -k1,1V -k2,2g -k3,3g {0} -o {0}".format(output_file))
    os.system("sort -k1,1V -k2,2g -k3,3g {0} -o {0}".format(ucscoutput_file))

    print("\n- Your text  output file is: \n\t- ENSEMBL: {0}\n\t- UCSC   : {1}".format(output_file, ucscoutput_file))
    # head /home/rad/users/gaurav/projects/pdacMetastasis/output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2.bed
    # ┌───┬──────────┬──────────┬────────────────────────────────────────┐
    # │ 1 │ 46083692 │ 46083693 │ DS08_5320_LivMet-1|1_46083693_T_G|8|no │
    # │ 1 │ 66632103 │ 66632104 │ DS08_5320_LivMet-1|1_66632104_C_T|4|no │
    # │ 1 │ 74099410 │ 74099411 │ DS08_5320_LivMet-1|1_74099411_G_C|1|no │
    # │ 1 │ 77514850 │ 77514851 │ DS08_5320_LivMet-1|1_77514851_A_C|3|no │
    # └───┴──────────┴──────────┴────────────────────────────────────────┘
    # ...


################ USER DEFINED FUNCTIONS ###################

def print_help():
    ''' Print system help '''
    print >> sys.stderr, "\n ----------------- HELP ------------------\n", parser.print_help(), "\n"

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python scripts/export_variant_input_to_bed.py -if=output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2.txt -od=/home/rad/users/gaurav/projects/pdacMetastasis/output/mPDAC_cherrypicked1
        -------------------------------------------------
        -------------------------------------------------
        '''))

    # Add arguments
    parser.add_argument("-if", metavar='--infile', help="*Input file"                , dest="inputVariantFile"     , type=str, required=True)
    parser.add_argument("-of", metavar='--otfile', help="*Interim file path and name", dest="interim_path_basename", type=str, required=True)

    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().inputVariantFile:
        logdir="{0}/interimFiles/logs".format(get_file_info(get_file_info(parser.parse_args().inputVariantFile)[0])[0])
        create_dir(logdir)
        logfile = "{0}/{1}_export_variant_input_to_bed.log".format(logdir, get_file_info(parser.parse_args().inputVariantFile)[1])
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
    inputVariantFile      = args.inputVariantFile
    interim_path_basename = args.interim_path_basename

    # Get output file name
    output_file     = "{0}.bed".format(interim_path_basename)
    ofile           = open(output_file, 'w')

    ucscoutput_file = "{0}.ucscbed".format(interim_path_basename)
    uofile          = open(ucscoutput_file, 'w')

    main()
