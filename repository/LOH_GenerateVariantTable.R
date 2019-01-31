#!/usr/bin/Rscript

##########################################################################################
##
## LOH_GenerateVariantTable.R
##
## Calculate datapoints needed for plotting LOH while transforming them to B-allele frequencies.
##
##########################################################################################

args = commandArgs(TRUE)

name=args[1] #used for naming in- and output files
genome_fasta=args[2] #location of the reference genome fasta
repository_dir=args[3] #location of repository

source(paste(repository_dir,"/LOH_Library.R",sep=""))

#read input files
tumor = read.table(paste(name,"/results/Mutect2/",name,".Tumor.Mutect2.Positions.txt", sep=""), header=T, sep="\t")
normal = read.table(paste(name,"/results/Mutect2/",name,".Normal.Mutect2.Positions.txt", sep=""), header=T, sep="\t")

#adjust column names
colnames(tumor) = c("Chrom", "Pos", "Ref", "Alt", "Frequency", "RefCount", "AltCount", "MapQ", "BaseQ", "Strand_Art_F", "Strand_Art_R")
colnames(normal) = c("Chrom", "Pos", "Ref", "Alt", "Frequency", "RefCount", "AltCount", "MapQ", "BaseQ", "Strand_Art_F", "Strand_Art_R")

#filter reads
tumor=LOH_FilterReads(tumor)
normal=LOH_FilterReads(normal)

#merge both files using UniquePos, reduce the input to heterozygous positions
variants = LOH_MergeVariants(tumor,normal)

variants = variants[variants[,"Normal_Freq"] <= 0.7 & variants[,"Normal_Freq"] >= 0.3,]

#define the dictionary, which defines whether an allele im AMBigous or UNAMBigous.
dicts = LOH_DefineDictionaries()

dict_for=dicts[["dict_for"]]
dict_rev=dicts[["dict_rev"]]
dict_unamb=dicts[["dict_unamb"]]

#determine which variants can be assigned A/B immediately (unamb), and which need further lookup (amb)
results = LOH_AssignStatus(variants,dict_unamb)

amb = results[["amb"]]
indel = results[["indel"]]
unamb = results[["unamb"]]

#determine which allele serves as A- and B-allele
LOH_GenerateLookupTable(name,amb)

#extract positions surrounding the variants from the reference genome. These are needed for the assignemnt of the A- and B-Allele.
LOH_ExtractSurroundingNucleotides(name,genome_fasta)

#Re-Import all nucleotides, which surround the requested variants
import = LOH_ImportSurroundingNucleotides(name)

#Using the information from the surrounding nucleotides, all variants in amb are assigned to A or B
amb = LOH_AssignAlleles(amb,import,dict_unamb)

#merge the final tables and export it for plotting
final_variants = rbind(amb,unamb,indel)

write.table(final_variants,file=paste(name,"/results/LOH/",name,".VariantsForLOH.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)