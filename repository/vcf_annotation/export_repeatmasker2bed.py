#!/usr/local/bin/python
# coding: utf-8
"""
***********************************************
- PROGRAM: export_repeatmasker2bed.py
***********************************************
"""

def main():
    print("- Reading and converting ucsc repeatmasker text file to bed file...")
    with open(input_file, 'r') as f:
        # column -t input/annotation/mm10/mm10_ucsc_repeatMasker.txt | head
        # #swScore  genoName              genoStart  genoEnd    strand  repName            repClass        repFamily       repStart  repEnd  repLeft
        # 239       chr1                  67108752   67108881   +       RLTR17B_Mm         LTR             ERVK            329       450     -352
        # 230       chr1                  134217651  134217732  -       BC1_Mm             scRNA           scRNA           -56       96      14
        # ...

        # Get the header
        hline = "#CHROM\tREP_GENO_START\tREP_GENO_END\tRepName|Class|Family\tswScore\tSTRAND"
        ofile.write("{0}\n".format(hline))

        for line in (l.strip() for l in f if not l.startswith('#')):
            # Split the line
            row  = line.split("\t")

            # Create the name filed with all the columns of the input file separated by pipe
            name = "|".join([row[5],row[6],row[7]])

            # Output line in bed file (Ensembl format: No chr string)
            bedl = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(row[1].replace('chr',''), row[2], row[3], name, row[0], row[4])

            # Write to the output line
            ofile.write("{0}\n".format(bedl))

    ofile.close()
    print("\n- Your text  output file is: {0}".format(output_file))

    # Remove the non-chromosomes contigs
    print("egrep -v 'Un|random|GL|JH|M' {0} > {1} && mv {1} {0}\n".format(output_file,"{0}.tmp".format(output_file)))
    os.system("egrep -v 'Un|random|GL|JH|M' {0} > {1} && mv {1} {0}\n".format(output_file,"{0}.tmp".format(output_file)))


    # head /home/rad/users/gaurav/projects/pdacMetastasis/input/annotation/mm10/mm10_repeatMasker.bed
    # ┌───────┬────────────────┬──────────────┬──────────────────────┬─────────┬────────┐
    # │ CHROM │ REP_GENO_START │ REP_GENO_END │ RepName|Class|Family │ swScore │ STRAND │
    # ├───────┼────────────────┼──────────────┼──────────────────────┼─────────┼────────┤
    # │     1 │       67108752 │     67108881 │ RLTR17B_Mm|LTR|ERVK  │     239 │ +      │
    # │     1 │      134217651 │    134217732 │ BC1_Mm|scRNA|scRNA   │     230 │ -      │
    # │     1 │        8386825 │      8389555 │ Lx2|LINE|L1          │    8310 │ -      │
    # └───────┴────────────────┴──────────────┴──────────────────────┴─────────┴────────┘
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
        - python scripts/export_repeatmasker2bed.py -if=input/annotation/mm10/mm10_ucsc_repeatMasker.txt
        -------------------------------------------------
        -------------------------------------------------
        '''))

    # Add arguments
    parser.add_argument("-if", metavar='--infile', help="input  file", dest="input_file", type=str, required=True)

    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().input_file:
        logdir="{0}/logs".format(get_file_info(get_file_info(parser.parse_args().input_file)[0])[0])
        create_dir(logdir)
        logfile = "{0}/{1}_export_repeatmasker2bed.log".format(logdir, get_file_info(parser.parse_args().input_file)[1])
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
    input_file = args.input_file

    # Get output file name
    output_file = "{0}.bed".format(get_file_info(input_file)[3]).replace('_ucsc','')
    ofile = open(output_file, 'w')

    main()
