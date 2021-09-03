#!/usr/local/bin/python
# coding: utf-8
"""
***********************************************
- PROGRAM: convert_vcf_to_annotation_input.py
***********************************************
"""


def main():

    # 1) Import original mergedvcf file
    # input_file = '/media/rad/HDD1/pdacMetastasis/input/mPDAC_wgs/DS4/results/Mutect2/DS4.Mutect2.txt'
    # ┌───────────────────┬───────────────────────┬─────────────────────────┬────────────────────┬──────┐
    # │    Old_Header     │      New_Header       │          SNV1           │        SNV2        │ SNV… │
    # ├───────────────────┼───────────────────────┼─────────────────────────┼────────────────────┼──────┤
    # │ CHROM             │ Chrom                 │ 1                       │ 1                  │ …    │
    # │ POS               │ GenomicPos            │ 37186919                │ 37186919           │ …    │
    # │ REF               │ Ref                   │ C                       │ C                  │ …    │
    # │ ALT               │ Alt                   │ T                       │ T                  │ …    │
    # │ GEN[Tumor].AF     │ TumorVAFMutect2       │ 0.402                   │ 0.402              │ …    │
    # │ GEN[Tumor].AD[0]  │ TumorRefCov           │ 69                      │ 69                 │ …    │
    # │ GEN[Tumor].AD[1]  │ TumorAltCov           │ 46                      │ 46                 │ …    │
    # │ GEN[Normal].AD[0] │ NormalRefCov          │ 120                     │ 120                │ …    │
    # │ GEN[Normal].AD[1] │ NormalAltCov          │ 0                       │ 0                  │ …    │
    # │ ANN[*].GENE       │ GeneName              │ Vwa3b                   │ Vwa3b              │ …    │
    # │ ANN[*].EFFECT     │ VariantLocationInGene │ downstream_gene_variant │ intragenic_variant │ …    │
    # │ ANN[*].IMPACT     │                       │ MODIFIER                │ MODIFIER           │ …    │
    # │ ANN[*].FEATUREID  │ GeneID                │ ENSMUST00000027289.13   │ ENSMUSG00000050122 │ …    │
    # │ ANN[*].HGVS_C     │ HGVSc                 │ c.*44154C>T             │ n.37186919C>T      │ …    │
    # │ ANN[*].HGVS_P     │ HGVSp                 │                         │                    │ …    │
    # └───────────────────┴───────────────────────┴─────────────────────────┴────────────────────┴──────┘

    print("- Reading and converting mergedvcf file to bed file...")
    mergedvcfDF = pd.read_csv(input_file, sep='\t', quoting=3)

    # Dropping extra columns
    print("\n- Dropping extra columns:\n\t- ANN[*].IMPACT")
    mergedvcfDF.drop(columns=['ANN[*].IMPACT'], inplace=True)

    # Rename existing columns
    print("\n- Renaming existing columns:")
    colNamesDict = {'CHROM':'Chrom', 'POS':'GenomicPos', 'REF':'Ref', 'ALT':'Alt', 'GEN[Tumor].AF':'TumorVAFMutect2', 'GEN[Tumor].AD[0]':'TumorRefCov', 'GEN[Tumor].AD[1]':'TumorAltCov', 'GEN[Normal].AD[0]':'NormalRefCov', 'GEN[Normal].AD[1]':'NormalAltCov', 'ANN[*].GENE':'GeneName', 'ANN[*].EFFECT':'VariantLocationInGene', 'ANN[*].FEATUREID':'GeneID', 'ANN[*].HGVS_C':'HGVSc', 'ANN[*].HGVS_P':'HGVSp'}
    for k, v in colNamesDict.items():
        print("\t- {0: <17}: {1}".format(k,v))
    mergedvcfDF.rename(columns=colNamesDict, errors="raise", inplace=True)

    # Create additional columns
    print("\n- Creating additional columns:\n\t- FileID\n\t- VariantID\n\t- SearchID\n\t- TumorTotalCov\n\t- TumorVAFCustom\n\t- Rescued")
    mergedvcfDF['FileID']         = mergedvcfDF.apply(lambda row: "{0}".format(get_file_info(input_file)[1].replace(".Mutect2","")), axis=1)
    mergedvcfDF['VariantID']      = mergedvcfDF.apply(lambda row: "{0}_{1:.0f}_{2}_{3}".format(row['Chrom'], row['GenomicPos'], row['Ref'], row['Alt']), axis=1)
    mergedvcfDF['SearchID']       = mergedvcfDF.apply(lambda row: "{0}:{1:.0f}".format(row['Chrom'], row['GenomicPos']), axis=1)
    mergedvcfDF['TumorTotalCov']  = mergedvcfDF['TumorRefCov'] + mergedvcfDF['TumorAltCov']
    mergedvcfDF['TumorVAFCustom'] = mergedvcfDF['TumorAltCov'] / mergedvcfDF['TumorTotalCov']
    mergedvcfDF['Rescued']        = "PPTonly"

    # Optional filtering
    mergedvcfDF = mergedvcfDF[mergedvcfDF['NormalAltCov'] == 0]

    # Reorder columns to fit the rescue format
    print("\n- Reordering columns to fit the rescue format")
    print("\t- FileID", "VariantID","SearchID", "Chrom", "GenomicPos", "Ref", "Alt", "TumorVAFMutect2", "TumorVAFCustom", "TumorTotalCov", "TumorRefCov", "TumorAltCov", "NormalRefCov", "NormalAltCov", "GeneName", "VariantLocationInGene", "GeneID", "HGVSc", "HGVSp", "Rescued")
    mergedvcfDF = mergedvcfDF[["FileID", "VariantID","SearchID", "Chrom", "GenomicPos", "Ref", "Alt", "TumorVAFMutect2", "TumorVAFCustom", "TumorTotalCov", "TumorRefCov", "TumorAltCov", "NormalRefCov", "NormalAltCov", "GeneName", "VariantLocationInGene", "GeneID", "HGVSc", "HGVSp", "Rescued"]]

    # Save the dataframe to output bed file (Ensembl format without chr)
    print("\n- Saving the results to the output file")
    mergedvcfDF.to_csv(output_file, sep="\t", index=False, float_format='%.3g')

    print("\n- Your text  output file is: \n\t- {0}".format(output_file))
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
    # │ totalCov       │ TumorTotalCov         │ 62                   │ 20                    │ …    │
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
        - python scripts/convert_vcf_to_annotation_input.py -if=/media/rad/HDD1/pdacMetastasis/input/mPDAC_wgs/DS4/results/Mutect2/DS4.Mutect2.txt -od=/home/rad/users/gaurav/projects/pdacMetastasis/output/mPDAC_wgs/txt
        -------------------------------------------------
        -------------------------------------------------
        '''))

    # Add arguments
    parser.add_argument("-if", metavar='--infile', help="*input file"     , dest="input_file", type=str, required=True)
    parser.add_argument("-od", metavar='--outdir', help="*Output directry", dest="output_dir", type=str, required=True)
    parser.add_argument("-sn", metavar='--smnam', help="*Sample name", dest="sampleName", type=str, required=True)

    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().input_file:
        logdir="{0}/interimFiles/logs".format(get_file_info(parser.parse_args().output_dir)[0])
        create_dir(logdir)
        logfile = "{0}/{1}_convert_vcf_to_annotation_input.log".format(logdir, get_file_info(parser.parse_args().input_file)[1])
        print(logdir)
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
  input_file      = args.input_file
  output_dir      = args.output_dir
  sampleName      = args.sampleName


  # Get output file name
  create_dir(output_dir)
  output_file     = "{0}/{1}_Mutect2.txt".format(output_dir, sampleName)
  ofile           = open(output_file, 'w')

  main()
