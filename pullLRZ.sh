# script to pull git repo on LRZ and set permissions automatically
git pull

# set permissions for self
chmod 754 pullLRZ.sh

# set permissions for main .sh files
chmod 754 MoCaSeq.sh MoCaSeq_LRZ_remap.sh MoCaSeq_LRZ_calling.sh ccc_MutectPostprocessingRunner.sh encrypt_LRZ_results.sh processResultsLRZ.sh

# set permissions for files in repository
chmod 754 repository/all_RunBubbleTree.R repository/SNV_SelectOutput.R repository/SNV_Signatures.R repository/CNV_PlotCNVKit.R repository/CNV_PlotHMMCopy.R repository/SV_MantaPostprocessing.sh repository/CNV_PlotCopywriter.R repository/LOH_Library.R repository/LOH_GenerateVariantTable.R repository/CNV_EstimateCoverage.R repository/CheckReferenceFiles.sh repository/CNV_RunCopywriter.R repository/Preparation_GenerateBWAIndex.sh repository/all_MoCaSeq_lcWGS.sh repository/vcf_annotation/get_usable_indels_for_annoatation.py repository/vcf_annotation/merge_annotation_files.py repository/vcf_annotation/export_variant_input_to_bed.py repository/vcf_annotation/export_repeatmasker2bed.py repository/vcf_annotation/parse_extended_mpileup.py repository/vcf_annotation/filter_variant_input.py repository/vcf_annotation/parse_mpileup.py repository/Preparation_GetExemplaryData.sh repository/LOH_MakePlots.R repository/CNV_RunHMMCopy.sh repository/SNV_SelectOutputSS.R

# set permissions for remap launcher files
chmod 754 launch/ccc_remap_wrapper.sh
chmod 754 launch/ccc_Mutect_Postprocessing.sh
# chmod 754 launch/remap/*.sh
chmod 754 launch/remap/batch0[1-4]-[A-Z][A-Z][A-Z][A-Z]*_[0-9][0-9][0-9][0-9]/*.sh

# set permissions for MoCaSeq launcher files
chmod 754 launch/mocaseq/batch0[1-4]-[A-Z][A-Z][A-Z][A-Z]*_[0-9][0-9][0-9][0-9]/*.sh

