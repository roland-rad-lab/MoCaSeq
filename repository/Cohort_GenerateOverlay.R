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

for (SummaryStat in c("Mean","Proportion"))
{
  for (Ylim in c(2,5))
  {
    RunOverlayAnalysis(Samples=Samples,Method="HMMCopy",species="Mouse",resolution=20000,SummaryStat=SummaryStat,AberrationCutoff=0.25,ChromomsomesToRemove="",Ylim=Ylim,format="pdf",Suffix=paste0(".",SummaryStat,".",Ylim))
  }
}

