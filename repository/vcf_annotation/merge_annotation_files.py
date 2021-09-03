#!/usr/local/bin/python -W ignore::DeprecationWarning
# coding: utf-8
"""
***********************************************
- PROGRAM: merge_annotation_files.py
***********************************************
"""

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
    mergedvcfDF = pd.read_csv(inputVariantFile, sep='\t', quoting=3); mergedvcfDF.name = 'mergedvcfDF'
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
    pileupDF = pd.read_csv(pileupFile, sep='\t', quoting=3); pileupDF.name = 'pileupDF'
    print("\t- Dropping redundant columns")
    print("\t\t- ['Ref', 'TotalCov']")
    pileupDF.drop(columns=['Ref', 'TotalCov'], inplace=True)

    # 3) Get the Segmental duplication annotated file
    print("\n- Reading the segmental duplication annotated file")

    # segDupFile = 'output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_segmentalDups.bed'
    # ┌───────┬──────────┬──────────┬────────────────────────────────────────┬────────────────────┐
    # │ Chrom │  Start   │GenomicPos│              VariantIDKey              │ SegmentalDupsCount │
    # ├───────┼──────────┼──────────┼────────────────────────────────────────┼────────────────────┤
    # │     1 │ 46083692 │ 46083693 │ DS08_5320_LivMet-1|1_46083693_T_G|8|no │                  1 │
    # │     1 │ 66632103 │ 66632104 │ DS08_5320_LivMet-1|1_66632104_C_T|4|no │                  0 │
    # │     1 │ 74099410 │ 74099411 │ DS08_5320_LivMet-1|1_74099411_G_C|1|no │                  0 │
    # │     1 │ 77514850 │ 77514851 │ DS08_5320_LivMet-1|1_77514851_A_C|3|no │                  0 │
    # │     1 │ 85098153 │ 85098154 │ DS08_5320_LivMet-1|1_85098154_G_T|3|no │                  2 │
    # └───────┴──────────┴──────────┴────────────────────────────────────────┴────────────────────┘
    segDupDF = pd.read_csv(segDupFile, sep='\t', quoting=3); segDupDF.name = 'segDupDF'
    print("\t- Dropping redundant columns")
    print("\t\t- ['Start', 'VariantIDKey']")
    segDupDF.drop(columns=['Start', 'VariantIDKey'], inplace=True)

    # 4) Get the Repeat masker annotated file
    print("\n- Reading the repeatmasker annotated file")
    # repMaskFile = 'output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_repeatmasker.txt'
    # ==================== ================== =======================================
    # Chrom                1                  1
    # Start                46083692           66632103
    # GenomicPos           46083693           66632104
    # VariantIDKey         1_46083693_T_G_8   1_66632104_C_T_4
    # RepChrom             .                  1
    # RepGenoStart         -1                 66632000
    # RepGenoEnd           -1                 66632176
    # RepNameClassFamily   .                  CT-rich|Low_complexity|Low_complexity
    # SwScore              -1                 885
    # SwStrand             .                  +
    # RepMaskFlag          0                  1
    # ==================== ================== =======================================
    repMaskDF = pd.read_csv(repMaskFile, sep='\t', quoting=3); repMaskDF.name = 'repMaskDF'
    print("\t- Dropping redundant columns")
    print("\t\t- ['Start', 'VariantIDKey']")
    repMaskDF.drop(columns=['Start', 'VariantIDKey'], inplace=True)

    # 5) Get the GC Percentage SNV base annotated file
    print("\n- Reading the GC Percentage SNV base annotated file")
    # gcpcSnvBaseFile = 'output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_5bp_GC_percentage.txt'
    # ========================= ================= ===================
    # GCpc5bpChrom              1                 1
    # GCpc5bpStart              7411766           7649717
    # GCpc5bpEnd                7411767           7649718
    # GCVariantBasePercentage   0                 80
    # Chrom                     1                 1
    # Start                     7411766           7649717
    # GenomicPos                7411767           7649718
    # VariantIDKey              1_7411767_A_G_0   1_7649718_CC_AT_2
    # ========================= ================= ===================
    gcpcSnvBaseDF = pd.read_csv(gcpcSnvBaseFile, sep='\t', usecols=['Chrom', 'GenomicPos', 'GCVariantBasePercentage'], quoting=3); gcpcSnvBaseDF.name = 'gcpcSnvBaseDF'

    # 6) Get the SNV reads mean GC percentage annotated file
    print("\n- Reading the SNV reads mean GC percentage annotated file")
    # gcpcSnvReadFile = 'output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_snvReads_Mean_GC_percentage.txt'
    # ====================== ================= =================== ================== ===================
    # Chrom                  1                 1                   1                  1
    # Start                  7411766           7649717             10554584           19228312
    # GenomicPos             7411767           7649718             10554585           19228313
    # VariantIDKey           1_7411767_A_G_0   1_7649718_CC_AT_2   1_10554585_C_G_0   1_19228313_AG_A_0
    # VariantReadsMeanGCPc   0.310191          0.437061            0.435377           0.50431
    # ====================== ================= =================== ================== ===================
    gcpcSnvReadDF = pd.read_csv(gcpcSnvReadFile, sep='\t', usecols=['Chrom', 'GenomicPos', 'VariantReadsMeanGCPc'], quoting=3); gcpcSnvReadDF.name = 'gcpcSnvReadDF'

    # 7) Get 10 base pairs around SNV annotated file
    print("\n- Reading the 10 base pairs around SNV annotated file")
    # tenBPAroundSNVFile = 'output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_10bp_around_SNV.txt'
    # ┌───────┬─────────┬─────────┬────────────┬─────────────────────────────┐
    # │ Chrom │  Start  │   End   │ GenomicPos │ Variant10bpFlankingSequence │
    # ├───────┼─────────┼─────────┼────────────┼─────────────────────────────┤
    # │     1 │ 7411756 │ 7411777 │    7411767 │ AGTCATGAGTATAAAGCCATG       │
    # │     1 │ 7649707 │ 7649728 │    7649718 │ atccacaaggccatggctgag       │
    # └───────┴─────────┴─────────┴────────────┴─────────────────────────────┘
    tenBPAroundSNVDF = pd.read_csv(tenBPAroundSNVFile, sep='\t', usecols=['Chrom', 'GenomicPos', 'Variant10bpFlankingSequence'], quoting=3); tenBPAroundSNVDF.name = 'tenBPAroundSNVDF'

    # 8) Get closest WES probe around SNV annotated file
    print("\n- Reading closest WES probe around SNV annotated file")
    # closestProbeFile = 'output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_closest_wes_probe.txt'
    # =========================== ================= ===================
    # Chrom                       1                 1
    # Start                       7411766           7649717
    # GenomicPos                  7411767           7649718
    # VariantIDKey                1_7411767_A_G_0   1_7649718_CC_AT_2
    # WESProbeChrom               1                 1
    # WESProbeStart               7177966           7177966
    # WESProbeEnd                 7178085           7178085
    # WESProbeAnnotation          ...               ...
    # WESProbeScore               1                 1
    # WESProbeStrand              -                 -
    # WESProbeDistanceToVariant   233682            471633
    # =========================== ================= ===================
    closestProbeDF = pd.read_csv(closestProbeFile, sep='\t', usecols=['Chrom', 'GenomicPos', 'WESProbeChrom', 'WESProbeStart', 'WESProbeEnd', 'WESProbeAnnotation', 'WESProbeScore', 'WESProbeStrand', 'WESProbeDistanceToVariant'], quoting=3); closestProbeDF.name = 'closestProbeDF'

    # 9) Get tumor indel file
    # head /home/rad/users/gaurav/projects/pdacMetastasis/output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_indels.txt
    # ┌───────┬────────────┬──────────────────────────────┬─────────────────────────┬────────────────────────────────────┐
    # │ Chrom │ GenomicPos │ TumorExactIndelSubsTypeCount │ TumorExactIndelSubsType │           InsertionBaseQ           │
    # ├───────┼────────────┼──────────────────────────────┼─────────────────────────┼────────────────────────────────────┤
    # │     1 │   88246255 │                            1 │ -7tgaaggt               │                                    │
    # │     1 │  191354613 │                            1 │ -5aagaa                 │                                    │
    # │     2 │  120720801 │                            1 │ +1g                     │ 33.29,0.86,31,33.00,33.00,34.00,34 │
    # └───────┴────────────┴──────────────────────────────┴─────────────────────────┴────────────────────────────────────┘
    print("\n- Reading tumor exact indel annotated file")
    tumorExactIndelDF = pd.read_csv(tumorExactIndelFile, sep='\t', quoting=3); tumorExactIndelDF.name = 'tumorExactIndelDF'

    # 10) Get tumor and normal extended  indel file
    # head /home/rad/users/gaurav/projects/pdacMetastasis/output/mPDAC_cherrypicked1/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_normal_parsed_normal_extendedpileup.txt
    # THIS IS TRANSPOSE VIEW OF FILE
    # ================================== ======================= ======================= =======================
    # Chrom                              1                       1                       2
    # GenomicPos                         88246255                191354613               120720801
    # SearchID                           1:88246255              1:191354613             2:120720801
    # NormalExtendedIndelTypeCount       0                       0                       0
    # NormalExtendedIndelTypeCountList   0|0|0|0|0|0|0|0|0|0|0   0|0|0|0|0|0|0|0|0|0|0   0|0|0|0|0|0|0|0|0|0|0
    # NormalExtendedIndelSubsType        ||||||||||              |||||g|||||             ||||||||||
    # NormalExtendedIndelSubsPos         …                       …                       …
    # NormalExtendedIndelSubsReadBases   …                       …                       …
    # ================================== ======================= ======================= =======================
    print("\n- Reading normal extended indel annotated file")
    tumorExtendedIndelDF = pd.read_csv(tumorExtendedIndelFile, sep='\t', quoting=3); tumorExtendedIndelDF.name = 'tumorExtendedIndelDF'
    print("\t- Dropping redundant columns")
    print("\t\t- ['SearchID']")
    tumorExtendedIndelDF.drop(columns=['SearchID'], inplace=True)

    print("\n- Reading normal extended indel annotated file")
    normalExtendedIndelDF = pd.read_csv(normalExtendedIndelFile, sep='\t', quoting=3); normalExtendedIndelDF.name = 'normalExtendedIndelDF'
    print("\t- Dropping redundant columns")
    print("\t\t- ['SearchID']")
    normalExtendedIndelDF.drop(columns=['SearchID'], inplace=True)

    ############################################################################################
    # Merge the DFs to the originalDF on ['Chrom', 'GenomicPos']
    # Source: https://stackoverflow.com/questions/54049458/merge-multiple-dataframes-by-multiple-id-columns-in-pandas
    print("\n- Merging the DFs to the originalDF on ['Chrom', 'GenomicPos']")
    print("\t- [mergedvcfDF, pileupDF, segDupDF, repMaskDF, gcpcSnvBaseDF, gcpcSnvReadDF, tenBPAroundSNVDF, closestProbeDF]")
    annotatedSnvDF = [mergedvcfDF, pileupDF, segDupDF, repMaskDF, gcpcSnvBaseDF, gcpcSnvReadDF, tenBPAroundSNVDF, closestProbeDF]

    #import pdb; pdb.set_trace()

    print("\n- Merge DFs")

    # NIKLAS QUICK FIX
    aDF = mergedvcfDF.drop_duplicates().copy(deep=True)
    bDF = pileupDF.drop_duplicates().copy(deep=True)
    cDF = segDupDF.drop_duplicates().copy(deep=True)
    dDF = repMaskDF.drop_duplicates().copy(deep=True)
    eDF = gcpcSnvBaseDF.drop_duplicates().copy(deep=True)
    fDF = gcpcSnvReadDF.drop_duplicates().copy(deep=True)
    gDF = tenBPAroundSNVDF.drop_duplicates().copy(deep=True)
    hDF = closestProbeDF.drop_duplicates().copy(deep=True)
    tumorExactIndelDF = tumorExactIndelDF.drop_duplicates().copy(deep=True)
    tumorExtendedIndelDF = tumorExtendedIndelDF.drop_duplicates().copy(deep=True)
    normalExtendedIndelDF = normalExtendedIndelDF.drop_duplicates().copy(deep=True)

    annotatedSnvDF = [aDF, bDF, cDF, dDF, eDF, fDF, gDF, hDF]
    for i in annotatedSnvDF: i['Chrom'] = i['Chrom'].astype(str)
    outer_merge = partial(pd.merge, how='left')
    annotatedSnvDF = reduce(outer_merge, annotatedSnvDF)

    # GAURAV's broken and not working shit code
    # for i in annotatedSnvDF: i['Chrom'] = i['Chrom'].astype(str)
    # outer_merge = partial(pd.merge, how='outer')
    # annotatedSnvDF = reduce(outer_merge, annotatedSnvDF)


    # Step 2) Merge indel annotation: tumorExactIndelDF, normalExtendedIndelDF
    # Convert the chromosome type to oject in case of X, Y and M chromosomes
    annotatedSnvDF['Chrom']        = annotatedSnvDF['Chrom'].astype(str)
    tumorExactIndelDF['Chrom']     = tumorExactIndelDF['Chrom'].astype(str)
    tumorExtendedIndelDF['Chrom']  = tumorExtendedIndelDF['Chrom'].astype(str)
    normalExtendedIndelDF['Chrom'] = normalExtendedIndelDF['Chrom'].astype(str)
    annotatedSnvDF = pd.merge(annotatedSnvDF, tumorExactIndelDF    , how='left', on=['Chrom', 'GenomicPos'])
    annotatedSnvDF = pd.merge(annotatedSnvDF, tumorExtendedIndelDF , how='left', on=['Chrom', 'GenomicPos'])
    annotatedSnvDF = pd.merge(annotatedSnvDF, normalExtendedIndelDF, how='left', on=['Chrom', 'GenomicPos'])


    #import pdb; pdb.set_trace()

    # Add the insertion basequality to the meanBaseQuality columns
    # annotatedSnvDF['BaseQMean'].fillna(annotatedSnvDF['BaseQInsertion'], inplace=True)
    # annotatedSnvDF.drop(columns=['BaseQInsertion'], inplace=True)
    try:
        insbqDF = annotatedSnvDF['BaseQInsertion'].str.split(",", expand=True)
        annotatedSnvDF['BaseQMean'  ].fillna(insbqDF[0], inplace=True)
        annotatedSnvDF['BaseQStd'   ].fillna(insbqDF[1], inplace=True)
        annotatedSnvDF['BaseQMin'   ].fillna(insbqDF[2], inplace=True)
        annotatedSnvDF['BaseQ25pc'  ].fillna(insbqDF[3], inplace=True)
        annotatedSnvDF['BaseQMedian'].fillna(insbqDF[4], inplace=True)
        annotatedSnvDF['BaseQ75pc'  ].fillna(insbqDF[5], inplace=True)
        annotatedSnvDF['BaseQMax'   ].fillna(insbqDF[6], inplace=True)
    except:
        pass
    annotatedSnvDF.drop(columns=['BaseQInsertion'], inplace=True)

    # Check if detected IndelType is the Indel called by Mutect2
    # Sometimes Mutect2 detects e.g. an insertion at a specific site in one sample, but by chance there is at exact the same site a deletion occuring in another sample. This deletion would then wrongly pass the indelfilter, as it would only have one indeltype at that site.
    # E.g. 5_73064610_T_TG in DS08_5320_PPT-1
    # annotatedSnvDF['TumorExactIndelSubsType'].str.contains("[\+-]")
    # This Problem can be fixed by classifying a Mutect2 call as an Indel if: The sign(len(ALT) - len(REF)) [Ex: TG - T ] but the TumorExactIndelSubsType ==1 )has the opposite sign with the same length (-1g) then mark them as fail (0) else everything in the IndelCalledByMutect2_pass1_fail0 column is (1)
    # # Get the first character of the string in the column
    # annotatedSnvDF['TumorExactIndelSubsType'].str[0]
    # # Map a dictionary to get the numeric sign for the +/-
    # di = {'+': 1, '-': -1}
    # annotatedSnvDF['TumorExactIndelSubsType'].str[0].map(di)
    # Check if the sign of diffence of length of ALT-REF is same or different
    annotatedSnvDF['IndelCalledByMutect2Pass1Fail0'] = np.where(np.logical_and(annotatedSnvDF['TumorExactIndelSubsTypeCount']==1, (annotatedSnvDF['TumorExactIndelSubsType'].str[0].map({'+': 1, '-': -1}) != np.sign(annotatedSnvDF['Alt'].str.len()-annotatedSnvDF['Ref'].str.len()))),0,1)

    # Add final indel annoation column (PASS (1) or FAIL (0))
    # https://stackoverflow.com/questions/31413286/how-to-create-new-column-in-a-df-based-on-multiple-conditions
    annotatedSnvDF['IndelQualityPass1Fail0'] = np.where(np.logical_and(annotatedSnvDF['IndelCalledByMutect2Pass1Fail0'] ==1, np.logical_and(annotatedSnvDF['TumorExactIndelSubsTypeCount']==1, np.logical_and(annotatedSnvDF['TumorExtendedIndelTypeCount']==1, annotatedSnvDF['NormalExtendedIndelTypeCount']==0))),1,0)

    # Get the size of input and output DF
    print("\n- Input and output DF size:\n\t- Input  DF: {0}\n\t- Output DF: {1}".format(mergedvcfDF.shape, annotatedSnvDF.shape))
    nrowsBefore = annotatedSnvDF.shape[0]

    # Dropping duplicate values
    print("\n- Dropping duplicate values")
    annotatedSnvDF.drop_duplicates(subset='VariantID', inplace=True)

    # Number of duplicates dropped
    print("\t- Number of duplicates dropped: {0}".format(nrowsBefore - annotatedSnvDF.shape[0]))

    # Get the size of input and output DF
    print("\n- Input and output DF size:\n\t- Input  DF: {0}\n\t- Output DF: {1}".format(mergedvcfDF.shape, annotatedSnvDF.shape))
    # x = set(mergedvcfDF['VariantID']).symmetric_difference(set(annotatedSnvDF['VariantID']))
    # print(x)

    # Get the additional columns into a dataframe
    additionalDF = annotatedSnvDF[['VariantID','SearchID','ReadBases','BaseQ','MapQ','ReadName','TumorExtendedIndelSubsReadBases','NormalExtendedIndelSubsReadBases']]

    # Drop additinal columns
    annotatedSnvDF.drop(columns  =['ReadBases','BaseQ','MapQ','ReadName','TumorExtendedIndelSubsReadBases','NormalExtendedIndelSubsReadBases'], inplace=True)

    # Replace word 'nan' to empty string
    annotatedSnvDF = annotatedSnvDF.replace('nan', '', regex=True)

    #import pdb; pdb.set_trace()

    # Save to text file
    print("\n- Saving the results to the output files:")
    print("\t- Saving to tab separated text file")

    output_dir = "{0}/finalAnnotation/".format(get_file_info(inputVariantFile)[0])
    create_dir(output_dir)

    # old wrong code from Gaurav
    #output_txtfile = "{0}/finalAnnotation/{1}_variant_annotation.txt".format(get_file_info(get_file_info(inputVariantFile)[0])[0],get_file_info(inputVariantFile)[1])
    output_txtfile = "{0}/finalAnnotation/{1}_variant_annotation.txt".format(get_file_info(inputVariantFile)[0],get_file_info(inputVariantFile)[1])
    annotatedSnvDF.to_csv(output_txtfile, sep="\t", index=False)

    # Save to summary excel
    print("\t- Saving to MS excel file")
    #output_excel   = "{0}/finalAnnotation/{1}_variant_annotation.xlsx".format(get_file_info(get_file_info(inputVariantFile)[0])[0],get_file_info(inputVariantFile)[1])
    output_excel   = "{0}/finalAnnotation/{1}_variant_annotation.xlsx".format(get_file_info(inputVariantFile)[0],get_file_info(inputVariantFile)[1])
    create_dir(get_file_info(output_excel)[0])
    annotatedSnvDF.to_excel(output_excel, sheet_name='snvAnnotation', index=False, float_format="%.2g")

    # Save additional columns to text file
    print("\t- Saving additional meta data columns to tab separated text file")
    #addout_txtfile = "{0}/interimFiles/{1}/{2}_variant_annotation.txt".format(get_file_info(get_file_info(inputVariantFile)[0])[0],sampleName, get_file_info(inputVariantFile)[1])
    addout_txtfile = "{0}/interimFiles/{1}/{2}_variant_annotation.txt".format(get_file_info(inputVariantFile)[0],sampleName, get_file_info(inputVariantFile)[1])
    additionalDF.to_csv(addout_txtfile, sep="\t", index=False, header=True)

    # Print the path to output files
    print("\n- Output files:")
    print("\t- Your text  output file is: {0}".format(output_txtfile))
    print("\t- Your excel output file is: {0}".format(output_excel))

################ USER DEFINED FUNCTIONS ###################
def remove_duplicates_from_DFs(df):
    ''' Remove rows that have duplicate entries in the Chrom and POS columns'''
    return df[~df.duplicated(['Chrom','GenomicPos'])]

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python scripts/merge_annotation_files.py -mf=output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2.txt -pf=output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_readQ_mappingQuality.txt -sf=output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_segmentalDups.bed -rf=output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_repeatmasker.txt -gb=output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_5bp_GC_percentage.txt -ab=output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_10bp_around_SNV.txt -gr=output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_snvReads_Mean_GC_percentage.txt -cp=output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_closest_wes_probe.txt -ti=output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_indels.txt -tx=output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_parsed_tumor_extendedpileup.txt -nx=output/AGRad_mPDAC/mergedvcf/txt/DS08_5320_LivMet-1-Mutect2_parsed_normal_extendedpileup.txt
        -------------------------------------------------
        -------------------------------------------------
        '''))

    # Add arguments
    parser.add_argument("-sn", metavar='--sample', help="Name of the sample. ex: DS4"     , dest="sampleName"             , type=str, required=True)
    parser.add_argument("-if", metavar='--vafile', help="input variant annotation file"   , dest="inputVariantFile"       , type=str, required=True)
    parser.add_argument("-pf", metavar='--plfile', help="input pileup file"               , dest="pileupFile"             , type=str, required=True)
    parser.add_argument("-sf", metavar='--sgfile', help="input segmental dup file"        , dest="segDupFile"             , type=str, required=True)
    parser.add_argument("-rf", metavar='--rpfile', help="input repeat masker file"        , dest="repMaskFile"            , type=str, required=True)
    parser.add_argument("-gb", metavar='--gbfile', help="input GC %% for SNV base file"   , dest="gcpcSnvBaseFile"        , type=str, required=True)
    parser.add_argument("-gr", metavar='--grfile', help="input Mean GC %% for SNV Reads"  , dest="gcpcSnvReadFile"        , type=str, required=True)
    parser.add_argument("-cp", metavar='--cpfile', help="input closest WES probe for SNV" , dest="closestProbeFile"       , type=str, required=True)
    parser.add_argument("-ab", metavar='--abfile', help="input 10 bp around SNV file"     , dest="tenBPAroundSNVFile"     , type=str, required=True)
    parser.add_argument("-ti", metavar='--tifile', help="input tumor exact indel file"    , dest="tumorExactIndelFile"    , type=str, required=True)
    parser.add_argument("-tx", metavar='--nxfile', help="input tumor extended indel file" , dest="tumorExtendedIndelFile" , type=str, required=True)
    parser.add_argument("-nx", metavar='--nxfile', help="input normal extended indel file", dest="normalExtendedIndelFile", type=str, required=True)

	# http://thomas-cokelaer.info/blog/2014/03/python-argparse-issues-with-the-help-argument-typeerror-o-format-a-number-is-required-not-dict/

    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().inputVariantFile:
        logdir="{0}/interimFiles/logs".format(get_file_info(get_file_info(parser.parse_args().inputVariantFile)[0])[0])
        create_dir(logdir)
        logfile = "{0}/{1}_merge_annotation_files.log".format(logdir, get_file_info(parser.parse_args().inputVariantFile)[1])
    else:
        logdir  = "{0}/logs".format(os.getcwd())
        create_dir(logdir)
        logfile = "{0}/{1}_merge_annotation_files.log".format(logdir,get_file_info(sys.argv[0])[1])

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
    sampleName              = args.sampleName
    inputVariantFile        = args.inputVariantFile
    pileupFile              = args.pileupFile
    segDupFile              = args.segDupFile
    repMaskFile             = args.repMaskFile
    gcpcSnvBaseFile         = args.gcpcSnvBaseFile
    gcpcSnvReadFile         = args.gcpcSnvReadFile
    tenBPAroundSNVFile      = args.tenBPAroundSNVFile
    closestProbeFile        = args.closestProbeFile
    tumorExactIndelFile     = args.tumorExactIndelFile
    tumorExtendedIndelFile  = args.tumorExtendedIndelFile
    normalExtendedIndelFile = args.normalExtendedIndelFile

    main()
