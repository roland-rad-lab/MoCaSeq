#!/usr/local/bin/python
# coding: utf-8
"""
***********************************************
- PROGRAM: get_usable_indels_for_annoatation.py
- CONTACT: Gaurav Jain(gaurav.jain@tum.de)
***********************************************
"""
print(__doc__)

def main():

    # 1) Import original mergedvcf file
    # inputVariantFile = 'output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2.txt'
    # ======================= ======================
    # FileID                  DS08_5320_LivMet-1
    # VariantID               1_46083693_T_G
    # SearchID                1:46083693
    # Chrom                   1
    # GenomicPos              46083693
    # Ref                     T
    # Alt                     G
    # TumorVAFMutect2         0.129032
    # TumorVAFCustom          0.129032258064516
    # TotalCov                62
    # TumorRefCov             54
    # TumorAltCov             8
    # NormalRefCov            97
    # NormalAltCov            0
    # GeneName                Dnah7b
    # VariantLocationInGene   intron_variant
    # GeneID                  ENSMUST00000069293.8
    # HGVSc                   c.274-12T>G
    # HGVSp
    # Rescued                 no
    # ======================= ======================
    print("\n- Reading the base and mapping quality annotated file ...")
    mergedvcfDF = pd.read_csv(inputVariantFile, sep='\t', quoting=3)
    # Drop rows with NaN or missing value in all columns
    mergedvcfDF = mergedvcfDF.dropna(how='all')
    mergedvcfDF['GenomicPos'] = mergedvcfDF['GenomicPos'].apply(np.int64)

    # Remove non-standard chromosome rows
    try:
        mergedvcfDF = mergedvcfDF[~mergedvcfDF['Chrom'].str.contains('Un|random|GL|JH|M|KQ|K')]
    except:
        pass

    # 2) Get the base and mapping quality annotated file
    # pileupFile = 'output/mPDAC_cherrypicked1/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_parsed_tumor_exact_pileup.txt'
    # ================ ======================= ======================================================================
    # Chrom            1                       1
    # GenomicPos       7649718                 46083693
    # Ref              c                       t
    # TotalCov         20                      62
    # ReadBases        ,,.,.$A,,,,.,,,,,,a,,   ...........G.........G......,G.GG........G......GG.,...,..,^].^].^].
    # BaseQ            BBBBC@B@BCBCBCBCB?CB    @@@???>?<??B?>???????B??????@B?CC????????C??????CC?@??>@?@@??>
    # MapQ             >]SR@+J?.E;];9].L4@V    UUUUUUUUUUUTUUUUUUUUUTUUUUSUUY]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
    # ReadName         ...                     ...
    # BaseQMean        30.5                    33.62
    # BaseQStd         0.71                    0.52
    # BaseQMin         30                      33
    # BaseQ25pc        30.25                   33
    # BaseQMedian      30.5                    34
    # BaseQ75pc        30.75                   34
    # BaseQMax         31                      34
    # MapQMean         14.5                    57.25
    # MapQStd          6.36                    4.1
    # MapQMin          10                      51
    # MapQ25pc         12.25                   54.75
    # MapQMedian       14.5                    60
    # MapQ75pc         16.75                   60
    # MapQMax          19                      60
    # ================ ======================= ======================================================================
    pileupDF = pd.read_csv(parsedpileupFile, sep='\t', quoting=3)
    print("\t- Dropping redundant columns")
    print("\t\t- ['Ref', 'TotalCov']")
    pileupDF.drop(columns=['Ref', 'TotalCov'], inplace=True)
    pileupDF['Chrom']    = pileupDF['Chrom'].astype(str)
    mergedvcfDF['Chrom'] = mergedvcfDF['Chrom'].astype(str)

    # Remove duplicate entires from the dataframes
    print("\n- Merge the two input DFs and filter the entires ...")
    print("\t- Remove duplicate entires from the dataframes ...")
    mergedvcfDF = mergedvcfDF[~mergedvcfDF.duplicated(['Chrom','GenomicPos'])]
    pileupDF    = pileupDF[~pileupDF.duplicated(['Chrom','GenomicPos'])]

    # Merge the two input DFs and filter the entires
    print("\t- Merge DFs")
    annDF       = [mergedvcfDF, pileupDF]
    outer_merge = partial(pd.merge, how='outer')
    annDF       = reduce(outer_merge, annDF)

    # Drop rows with NaN or missing value in all columns
    print("\n- Drop rows with NaN or missing value in all columns")
    annDF = annDF.dropna(how='all')
    annDF['GenomicPos'] = annDF['GenomicPos'].apply(np.int64)
    # Remove non-standard chromosome rows
    print("\n- Remove non-standard chromosome rows (Un|random|GL|JH|M|KQ|K)")
    annDF = annDF[~annDF['Chrom'].str.contains('Un|random|GL|JH|M|KQ|K')]

    # Get the count of different INDEL types
    # itcln,_itln,ill,readBases = get_clean_readbase_indel_counts(row[4])
    print("\n- Get the count of different INDEL types")
    annDF['TumorExactIndelSubsTypeCount'], annDF['TumorExactIndelSubsType'], annDF['TumorExactIndelLocation'], _readBases = zip(*annDF['ReadBases'].apply(get_clean_readbase_indel_counts))

    # debugPos = '1:74273721'
    # print(annDF[annDF['SearchID']==debugPos].values.tolist())
    # pd.set_option('display.expand_frame_repr', False)
    # print(annDF['TumorExactIndelSubsTypeCount'][annDF['GenomicPos']==debugPos])

    # Filter out everthing except an indel
    # Variant positions: SNVs, DNVs, INDELs
    print("\n- Filter out all variant positions other than the INDELs")
    print("\t- Number of variant positions before filtering: {}".format(annDF.shape[0]))
    annDF  = annDF[annDF.apply(lambda row: len(row['Ref']) - len(row['Alt']), axis=1) != 0]
    print("\t- Number of variant positions that are INDELs : {}".format(annDF.shape[0]))


    # Filter out INDEs with 0 types
    print("\n-  Filter out INDEs with 0 types")
    # print(annDF.loc[38])
    annDF = annDF[annDF['TumorExactIndelSubsTypeCount'] != 0]
    annDF = annDF[annDF['TumorExactIndelSubsTypeCount'].notnull()]
    print("\t- Total variant positions that are actual INDELs : {}".format(annDF.shape[0]))
    # print(annDF)

    if(annDF.shape[0]):
        # Get the base quality of insersions
        # Get the read names for the reads containing the insertions (filter out reads with only deletions)
        insDF = annDF[(annDF['TumorExactIndelSubsType'].str.contains("[\+]")) & (~annDF['TumorExactIndelSubsType'].str.contains("[\-]"))]
        print("\t- Total variant positions that are INSERTIONS : {}".format(insDF.shape[0]))
        # print(insDF)
        insDF = get_insertion_baseQuality(subsetBamFile, insDF)

    # Add the insertion base quality column in annDF
    # annDF.loc[annDF['country'].isnull(),'country'] = insDF['issuer']
    try:
        annDF['BaseQInsertion'] = insDF['BaseQInsertion']
    except:
        annDF['BaseQInsertion'] = np.nan

    # Save the selected columns in the dataframe
    print("\n- Save the selected columns in the dataframe")
    annDF[['Chrom','GenomicPos', 'TumorExactIndelSubsTypeCount', 'TumorExactIndelSubsType', 'BaseQInsertion']].to_csv(output_file, sep="\t", index=False, float_format='%.f')

    # Create an empty output DF and add the relevant columns to it
    bedDF = pd.DataFrame(columns=['Chrom','start','end','chr_pos_ref_alt'])
    if(annDF.shape[0]):
        bedDF['Chrom'] = annDF['Chrom']
        bedDF['start'] = annDF['GenomicPos'] - 1
        bedDF['end']   = annDF['GenomicPos']
        bedDF['chr_pos_ref_alt']  = annDF.apply(lambda row: "{0}_{1:.0f}_{2}_{3}_{4:.0f}".format(row['Chrom'], row['GenomicPos'], row['Ref'], row['Alt'], row['TumorAltCov']), axis=1)

    # Save the dataframe to output bed file (Ensembl format without chr)
    print("\n- Saving the results to the output files")
    bedDF.to_csv(indel_file, sep="\t", index=False, header=False,  float_format='%.f')

    # Save the dataframe to ucscoutput bed file (UCSC bed format with chr)
    bedDF['Chrom'] = 'chr' + bedDF['Chrom'].astype(str)
    bedDF.to_csv(ucindel_file, sep="\t", index=False, header=False,  float_format='%.f')

    print("\n- Your text  output file is: \n\t- INDEL  : {0}\n\t- ENSEMBL: {1}\n\t- UCSC   : {2}".format(output_file, indel_file, ucindel_file))

################ USER DEFINED FUNCTIONS ###################

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
    subsList       = []
    indelType      = ''  # + or -
    indelLoc       = 0
    indelLocList   = []
    try:
        for i, x in enumerate(bases):
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
        # indelList = get_unique_list(list(chain.from_iterable(indelList))[1:])
        indelList = list(chain.from_iterable(indelList))[1:]

        # debugBases = '........-10TGTGTGTGTA..-10TGTGTGTGTA......'
        # debugBases = ',$..,..,..,..,.,.,.,,.,...,..,,...,-4tgtg,-10tgtgtgtgtg...-14TGTGTGTGTGTGTG,-2tg....,.,..*^I.^I,^]g'.upper()
        # if bases.upper() == debugBases:
        #     print("IndelList: {}".format(indelList))

        # This site has a 10 bp long Indel, but at the TumorExactIndelSubsType it is wrongly called as two Indels: -1|0tgtgtgtgta.
        # Correct would be: -10tgtgtgtgta
        # Location: 9_42008058_GTGTGTGTGTA_G
        # Merge the indexes correctly
        # input : ['+1a','-1', '0tgtgtgtgta','+2t', '-3','0ss']
        # output: ['+1a', '-10tgtgtgtgta', '+2t', '-30ss']
        idx    = []
        idxLoc = 0
        for i in indelList:
            # if bases.upper() == debugBases: print("{0}: {1}".format(i, i.islower()))
            if not i.islower():
                # idx.append(indelList.index(i))
                idx.append(idxLoc)
            idxLoc+=1

        # if bases == debugBases: print("idx:{}".format(idx))
        xl = indelList.copy()

        # Merge the elements
        for j in idx:
            # if bases == debugBases: print(xl[j:j+2])
            xl[j:j+1] = [''.join(xl[j:j+2])]
            # if bases == debugBases: print("after idx:{}".format(xl))

        # if bases == debugBases: print("after idx:{}".format(xl))

        # Remove the non-merged entries
        for i in sorted(idx, reverse=True):
            del xl[i+1]

        indelList = get_unique_list(xl.copy())
        # if bases == debugBases:
        #     print("final idx:{}".format(xl))
        #     try:
        #         sys.exit() # this always raises SystemExit
        #     except SystemExit:
        #         print("sys.exit() worked as expected")
        #         print(sys.exc_info()[0])
        #     except:
        #         print("Something went horribly wrong") # some other exception got raised

        # Get substitutions
        # 16_30547590_GCACA_G
        # bases=',$,,,,,,,,g,,,,.,,,,-2ca,-2ca.,,,,t'
        # subsList = ','.join(map(str, re.findall('[acgt]', callString))) # ['g', 't']
        subsList = get_unique_list(re.findall('[acgt]', callString)) # ['g', 't']

        # Add the substitutions along with the types of indels
        indelList.extend(subsList)

        # Get the length of the list to obtain number of unique indels
        # return(len(indelList), callString)
        return(int(len(indelList)), "|".join(indelList), indelLocList, callString)
    except:
        return(np.nan, np.nan, np.nan, np.nan)

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

def get_insertion_baseQuality(bamfile, insDF):
    '''
    Get base quality information for the insertions
    # https://sites.google.com/site/bioinformaticsremarks/bioinfo/sam-bam-format/what-is-a-cigar

    # 78       8_117024242_T_TTG         9                                          @@@?@@@@@  HS10_08370:1:1108:20102:5820#2,HS10_08370:1:21...
    # insDF[['ReadName']].loc[78].values.tolist()[0].split(",")
    # insDF[['ID', 'totalCov', 'BaseQuality', 'ReadName','TumorExactIndelSubsType','TumorExactIndelLocation']].loc[78]
        # ID8_117024242_T_TTG
        # totalCov9
        # BaseQuality@@@?@@@@@
        # ReadName                     HS10_08370:1:1108:20102:5820#2,HS10_08370:1:21...
        # TumorExactIndelSubsType+2tg
        # TumorExactIndelLocation                                            [5, 6, 8]

    # Get the reads containing indels
    # List:
        # In [57]: insDF[['ReadName']].loc[78].tolist()[0].split(",")
        # Out[57]:
        # ['HS10_08370:1:1108:20102:5820#2',
        #  'HS10_08370:1:2107:14360:11805#2',
        #  'HS10_08370:1:2110:6787:82707#2',
        #  'HS10_08370:1:2205:10169:5145#2',
        #  'HS10_08370:1:2215:20293:56551#2',
        #  'HS10_08370:1:2114:14943:69620#2',
        #  'HS10_08370:1:2302:21263:90486#2',
        #  'HS10_08370:1:1202:8940:65449#2',
        #  'HS10_08370:1:2106:7014:84227#2']
    #Indexes:
        # In [61]: insDF[['TumorExactIndelLocation']].loc[78].tolist()[0]
        # Out[61]: [5, 6, 8]

    # pd.Series(insDF[['ReadName']].loc[78].tolist()[0].split(","))[insDF[['TumorExactIndelLocation']].loc[78].tolist()[0]]
        # Out[67]:
        # 5    HS10_08370:1:2114:14943:69620#2
        # 6    HS10_08370:1:2302:21263:90486#2
        # 8     HS10_08370:1:2106:7014:84227#2
        # dtype: object

    # egrep 'HS10_08370:1:2114:14943:69620#2' /home/rad/users/gaurav/projects/pdacMetastasis/output/mPDAC_cherrypicked1/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2.sam
    # HS10_08370:1:2114:14943:69620#2	83	8	117024204	60	4S39M2I30M	=	117023986	-287	CATTCTTTCCCCCAACCCTGCTTCTAACCCAGTGACTAACCATTGTCTGGCCATTTCCCTCTGGACTCTAAGTTT	C@@;D@@;BBBBC@@BBBBBB@@C>@@BBCBABC?B>@@BC?@AA@CABBBC??@@BBB@CAAC@?@C>@CAA@?	SA:Z:JH590470.1,156859,+,75M,60,0;	MC:Z:67M7S	MD:Z:7G12G42T5	PG:Z:MarkDuplicates	RG:Z:1	NM:i:5	AS:i:46	XS:i:0	pa:f:0.613

    '''
    # Create a copy of original DF so that we avoid accidental errors by modifying original DF
    # https://www.dataquest.io/blog/settingwithcopywarning/
    insertionDF = insDF.copy()
    for idx in insertionDF.index.values:
        # Get the read names list for each insertion site
        # insReadNamesList = pd.Series(insertionDF[['ReadName']].loc[idx].tolist()[0].split(","))[insertionDF[['TumorExactIndelLocation']].loc[idx].tolist()[0]].tolist()
        insReadNamesList = pd.Series(insertionDF[['ReadName']].loc[idx].tolist()[0].split(",")).tolist()
        baseQualityList  = list()
        # print(insReadNamesList, baseQualityList)
        # print("For loop starts here...\n\n")

        # Get the bam file handle and loop through the file
        pysamhandle      = pysam.AlignmentFile(bamfile, 'rb')
        for read in pysamhandle:
            # Check if read name is present in the list
            if read.qname in insReadNamesList:
                cigarLine    = read.cigar
                # print(cigarLine)
                # To get the relative location and length of the insertion in the query string
                relLocation  = 0
                insertLength = 0
                isInsertion  = 0
                # For insertions get the mean basequality for each read and append it to the baseQualityList
                for cigarType, cigarLength in cigarLine:
                    # cigartuples: the cigar alignment. The alignment is returned as a list of tuples of (operation, length). If the alignment is not present, None is returned.
                    # The operations are:
                    # ┌───┬────────────────┬───┐
                    # │ M │ BAM_CMATCH     │ 0 │
                    # │ I │ BAM_CINS       │ 1 │ <----
                    # │ D │ BAM_CDEL       │ 2 │
                    # │ N │ BAM_CREF_SKIP  │ 3 │
                    # │ S │ BAM_CSOFT_CLIP │ 4 │
                    # │ H │ BAM_CHARD_CLIP │ 5 │
                    # │ P │ BAM_CPAD       │ 6 │
                    # │ = │ BAM_CEQUAL     │ 7 │
                    # │ X │ BAM_CDIFF      │ 8 │
                    # │ B │ BAM_CBACK      │ 9 │
                    # └───┴────────────────┴───┘
                    if cigarType == 1:
                        insertLength = cigarLength
                        isInsertion  = 1
                        break
                    else:
                        relLocation += cigarLength
                if isInsertion:
                    # print(cigarLine, relLocation, insertLength, read.query_sequence[relLocation:relLocation+insertLength], np.nanmean(read.query_qualities[relLocation:relLocation+insertLength]))
                    # print(baseQualityList)
                    insReadBaseQuality = np.mean(read.query_qualities[relLocation:relLocation+insertLength])
                    # print(read.qname, insReadBaseQuality)
                    baseQualityList.append(insReadBaseQuality)
                    insMappingQualityDF     = pd.DataFrame(baseQualityList).describe().T.values.tolist()[0]
                    insMappingQualityscores = "{0:1.2f},{1:1.2f},{2:1.0f},{3:1.2f},{4:1.2f},{5:1.2f},{6:1.0f}".format(insMappingQualityDF[1], insMappingQualityDF[2], insMappingQualityDF[3], insMappingQualityDF[4], insMappingQualityDF[5], insMappingQualityDF[6], insMappingQualityDF[7])

                    # Update the base quality scores with the new scores
                    # insertionDF.at[idx, 'BaseQInsertionDF'] = insMappingQualityscores
                    # print("*****************")
                    # print(insertionDF.shape)
                    # print("*****************")
                    insertionDF.loc[idx, 'BaseQInsertion'] = insMappingQualityscores

        # Set the insertion basequality value in a new column

    return insertionDF

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
        - python scripts/get_usable_indels_for_annoatation.py -mf=output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2.txt -pf=/home/rad/users/gaurav/projects/pdacMetastasis/output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_readQ_mappingQuality.txt -od=
        -------------------------------------------------
        CONTACT:
            Gaurav Jain
            gaurav.jain@tum.de
        -------------------------------------------------
        '''))

    # Add arguments
    parser.add_argument("-if", metavar='--vafile', help="*Input variant annotation file"   , dest="inputVariantFile"   , type=str, required=True)
    parser.add_argument("-pf", metavar='--plfile', help="*Input parsed pileup file", dest="parsedpileupFile", type=str, required=True)
    parser.add_argument("-bm", metavar='--bmfile', help="*Input indel subset bam file", dest="subsetBamFile", type=str, required=True)
    parser.add_argument("-od", metavar='--outdir', help="*Output directory", dest="output_dir", type=str, required=True)

    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().inputVariantFile:
        logdir="{0}/interimFiles/logs".format(get_file_info(get_file_info(parser.parse_args().inputVariantFile)[0])[0])
        create_dir(logdir)
        logfile = "{0}/{1}_get_usable_indels_for_annoatation.log".format(logdir, get_file_info(parser.parse_args().inputVariantFile)[1])
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
    from functools import partial, reduce

    # for looping files in a dir
    import glob

    # user defined modules
    from gjainPyLib import *      # import all the functions from the Gaurav`s python library

    ### for color scale
    from  matplotlib import colors
    from itertools import cycle, islice # barplot colors

    # For bam file
    import pysam

    ################ USER CONFIGURATION ###################
    np.set_printoptions(precision=6)
    #######################################################

    # Get input options
    args = check_options()

    # Store the variables
    parsedpileupFile = args.parsedpileupFile
    inputVariantFile = args.inputVariantFile
    subsetBamFile    = args.subsetBamFile
    output_dir       = args.output_dir

    # Get output file name
    output_file  = "{0}/{1}_indels.txt".format(output_dir    , get_file_info(inputVariantFile)[1])
    indel_file   = "{0}/{1}_indels.bed".format(output_dir    , get_file_info(inputVariantFile)[1])
    ucindel_file = "{0}/{1}_indels.ucscbed".format(output_dir, get_file_info(inputVariantFile)[1])

    main()
