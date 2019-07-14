#!/usr/bin/Rscript

##########################################################################################
##
## Cohort_GenerateOverlay.R
##
## Generates Overlay for all samples in the current directory.
##
##########################################################################################

repository_dir=args[1] #location of repository

source(paste(repository_dir,"/Cohort_GenerateOverlayLibrary.R",sep=""))     

Samples = list.dirs(full.name=F,recursive=F)
#Samples=Samples[-grep("multiqc_data",Samples)]
#Samples=Samples[-grep("REFERENCES",Samples)]
#Samples=Samples[-grep("1N7T85",Samples)]
#Samples=Samples[-grep("6P4EWP",Samples)]
#Samples=Samples[-grep("J837KK_Exo",Samples)]
#Samples=Samples[-grep("CohortDatabase",Samples)]
#Samples=Samples[-grep("CNVplot_combined_by_chr",Samples)]
#Samples=Samples[-grep("1_segment_annotated",Samples)]
#Samples=Samples[-grep("misc",Samples)]

for (SummaryStat in c("Mean","Proportion"))
{
  for (Ylim in c(2,5))
  {
    RunOverlayAnalysis(Samples=Samples,Method="HMMCopy",species="Human",resolution=10000,SummaryStat=SummaryStat,AberrationCutoff=0.25,ChromomsomesToRemove="",Ylim=Ylim,format="pdf",Suffix=paste0(".",SummaryStat,".",Ylim))
  }
}