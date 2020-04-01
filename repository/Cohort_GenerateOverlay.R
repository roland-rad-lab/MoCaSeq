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
Method="HMMCopy"
species="Mouse"
Paths=rep("",length(Samples))
Save_path=""

for (SummaryStat in c("Mean","Proportion"))
{
  for (Ylim in c(2,5))
  {
    RunOverlayAnalysis(Samples=Samples,Paths=Paths,Method=Method,species=species,resolution=20000,SummaryStat=SummaryStat,AberrationCutoff=0.25,ChromomsomesToRemove="",Ylim=Ylim,format="pdf",Suffix=paste0(".",SummaryStat,".",Ylim),Save_path=Save_path)
  }
}

