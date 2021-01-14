#!/usr/local/bin/python
# coding: utf-8
"""
***********************************************
- PROGRAM: filter_variant_input.py
- CONTACT: Gaurav Jain(gaurav.jain@tum.de)
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
    print("- Reading variant input file to datatable...")
    # mergedvcfDF = pd.read_csv(inputVariantFile, sep='\t', quoting=3)
    # Import data into the datatable
    start_time = time.time()
    variantDT = dt.fread(inputVariantFile, sep="\t", header=True, nthreads=16, quotechar="'",)
    print("%s seconds" % (time.time() - start_time))

    # Filter all rows with TumorVAFCustom == 0
    print("\t- Original variantDT shape: {0}".format(variantDT.shape))
    variantDT = variantDT[dt.f.TumorVAFCustom > 0,:]
    print("\t- Filtered variantDT shape: {0}".format(variantDT.shape))

    # Save the variantDT to output file
    variantDF = variantDT.to_pandas()
    variantDF.to_csv(output_file, sep='\t', index = False, float_format='%.4g')

    print("\n- Your filtered output file is: \n\t- {0}".format(output_file))
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
        - python scripts/filter_variant_input.py -if=/media/rad/HDD1/pdacMetastasis/output/mPDAC_wes/annInput/DS01_2259_LNMet-1_Mutect2.txt -of=/media/rad/HDD1/pdacMetastasis/output/mPDAC_wes/interimFiles/DS01_2259_LNMet-1/DS01_2259_LNMet-1_Mutect2
        -------------------------------------------------
        CONTACT:
            Gaurav Jain
            gaurav.jain@tum.de
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
    import datatable as dt
    from datatable import *
    import scanpy as sc
    import time

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
    output_file     = "{0}_vafFiltered.{1}".format(interim_path_basename, get_file_info(inputVariantFile)[2])
    create_dir(get_file_info(output_file)[0])
    ofile           = open(output_file, 'w')

    main()
