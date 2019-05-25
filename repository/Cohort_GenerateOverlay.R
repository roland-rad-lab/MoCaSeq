#!/usr/bin/Rscript

##########################################################################################
##
## Cohort_GenerateOverlay.R
##
## Generates Overlay for all samples in the current directory.
##
##########################################################################################

BinGenome = function(species,Method="Copywriter",resolution=10000,ChromomsomesToRemove=c("X","Y"))
              {
              Method <<- Method
              resolution <<- resolution
              ChromomsomesToRemove <<- ChromomsomesToRemove
              if(species=="Mouse")
                 {
                 ChromLength = c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,122082543,
                                 120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566,171031299,91744698)
                 names(ChromLength) = c(1:19,"X","Y")
                 ChromLength = ChromLength[!names(ChromLength) %in% ChromomsomesToRemove]
                 species <<- species
                 ChromLength <<- ChromLength                
                 }
               if(species=="Human")
                 {
                 ChromLength = c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,
                                 114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415)
                 names(ChromLength) = c(1:22,"X","Y")
                 ChromLength = ChromLength[!names(ChromLength) %in% ChromomsomesToRemove]
                 species <<- species
                 ChromLength <<- ChromLength                   
                 }                   
               BinnedGenome = data.frame(chr=NULL,start=NULL,end=NULL,strand=NULL,id=NULL)
               for ( chromosome in names(ChromLength))
                    {
                    Bins = seq(1,ChromLength[chromosome]-resolution,resolution)
                    BinsEnd = seq(resolution,ChromLength[chromosome],resolution)
                    BinnedChromosome = data.frame(chr=chromosome,start=Bins,end=BinsEnd)
                    BinnedGenome = rbind(BinnedGenome,BinnedChromosome)
                    }
              UnifiedPosition = paste(BinnedGenome$chr,BinnedGenome$start,BinnedGenome$end,sep="_")
              BinnedGenome$UnfiedPosition = UnifiedPosition
              BinnedGenome <<- BinnedGenome
              }




SetOverlayMat = function(Samples)
                      {
                      OverlayMat = matrix(0,nrow=nrow(BinnedGenome),ncol=3+(length(Samples)*2))
                      rownames(OverlayMat) = rownames(BinnedGenome)
                      colnames(OverlayMat) = c("Chromosome","Start","Stop",paste(Samples,"Gain",sep=""),paste(Samples,"Del",sep=""))
                      OverlayMat = as.data.frame(OverlayMat)
                      OverlayMat[,"Chromosome"] = BinnedGenome$chr
                      OverlayMat[,"Start"] = BinnedGenome$start
                      OverlayMat[,"Stop"] = BinnedGenome$end
                      OverlayMat <<- OverlayMat
                      }


####### Helper function  ###########
## Rename Copywriter Output       ##
## Rename Segment columns         ##
####################################

ReformatSegmentData = function(segmentdata)
                          {
                          method <<- Method
                          if(method=="Copywriter" & species=="Mouse")
                              {
                              segmentdata[segmentdata$Chrom==20,"Chrom"] = "X"
                              segmentdata[segmentdata$Chrom=="21","Chrom"] = "Y"
                              }
                          if(method=="Copywriter" & species=="Human")
                              {
                              segmentdata[segmentdata$Chrom==23,"Chrom"] = "X"
                              segmentdata[segmentdata$Chrom=="24","Chrom"] = "Y"                              
                              }
                          colnames(segmentdata) = c("chr","start","end","Mean")
                          UnifiedPosition = paste(as.factor(segmentdata$chr),segmentdata$start,floor(segmentdata$end),sep="_")
                          segmentdata$UnifiedPosition = UnifiedPosition
                          return(segmentdata)
                          }


LoadSegmentData = function(Sample)
                     {
                     MainPath = getwd()
                     if(Method=="Copywriter")
                         {
                         PathSuffix = "/results/Copywriter/"
                         SampleSuffix = ".Copywriter.segments.Mode.txt"
                         }
                     if(Method=="HMMCopy")
                        {
                        PathSuffix = "/results/HMMCopy/"
                        SampleSuffix = paste(".HMMCopy.",resolution,".segments.txt",sep="")
                        }    
                     file = paste(MainPath,"/",Sample,PathSuffix,Sample,SampleSuffix,sep="")
                     Segments = read.table(file,header=T,sep="\t")
                     return(Segments)
                     }




fillOverlayMat = function(Samples,AberrationCutoff,SummaryStat=SummaryStat)
                   {
                   SummaryStat <<- SummaryStat
                   for ( Sample in Samples)
                       {
                       Segments = LoadSegmentData(Sample)
                       Segments = ReformatSegmentData(Segments) 
                       Segments = Segments[abs(Segments$Mean)>=AberrationCutoff,]
                       Segments = Segments[!Segments$chr %in% ChromomsomesToRemove,]
                       for (segment in Segments$UnifiedPosition)
                            {
                            ActualPosition = Segments[Segments$UnifiedPosition== segment,]
                            Chromosome = ActualPosition$chr
                            Start = ActualPosition$start
                            Stop  = ActualPosition$end
                            if(ActualPosition$Mean>0)
                                {
                                OverlayMat[OverlayMat[,"Chromosome"]==Chromosome & OverlayMat[,"Start"]>=Start & OverlayMat[,"Stop"]<=Stop,paste(Sample,"Gain",sep="")] = ActualPosition$Mean
                                }
                            if(ActualPosition$Mean<0)
                                {
                                OverlayMat[OverlayMat[,"Chromosome"]==Chromosome & OverlayMat[,"Start"]>=Start & OverlayMat[,"Stop"]<=Stop,paste(Sample,"Del",sep="")] = ActualPosition$Mean
                                }
                            }
                       }
                   ColumnsGainInfo = grep("Gain",colnames(OverlayMat),value=T)
                   ColumnsDelInfo = grep("Del",colnames(OverlayMat),value=T)
                   if(SummaryStat=="Mean")
                       {
                       OverlayMat$Gain = apply(OverlayMat[,ColumnsGainInfo],1,mean)
                       OverlayMat$Del = apply(OverlayMat[,ColumnsDelInfo],1,mean)
                       }
                   if(SummaryStat=="Proportion")
                       {
                       AnnotationInfo = OverlayMat[,c("Chromosome","Start","Stop")]
                       CNVInfo = OverlayMat[,c(ColumnsGainInfo,ColumnsDelInfo)]
                       CNVInfo[CNVInfo!=0] = 1
                       OverlayMat = cbind(AnnotationInfo,CNVInfo)
                       OverlayMat$Gain = apply(OverlayMat[,ColumnsGainInfo],1,sum)/length(ColumnsGainInfo)
                       OverlayMat$Del = -1* apply(OverlayMat[,ColumnsDelInfo],1,sum)/length(ColumnsGainInfo)
                       }    
                   OverlayMat <<- OverlayMat
                   }    


ReformatOverlayMat = function()
                      {
                      if("X" %in% names(ChromLength) & "Y" %in% names(ChromLength))
                           {
                           Chromosomes = names(ChromLength)[!names(ChromLength) %in% c("X","Y")]
                           Chromosomes = c(sort(as.numeric(Chromosomes)),"X","Y")
                           }
                      if("X" %in% names(ChromLength) & !("Y" %in% names(ChromLength)))
                           {
                           Chromosomes = names(ChromLength)[names(ChromLength)!="X"]
                           Chromosomes = c(sort(as.numeric(Chromosomes)),"X")                           
                           }
                      if(!("X" %in% names(ChromLength)) & ("Y" %in% names(ChromLength)))
                           {
                           Chromosomes = names(ChromLength)[names(ChromLength)!="Y"]
                           Chromosomes = c(sort(as.numeric(Chromosomes)),"Y")
                           }                           
                      if(!("X" %in% names(ChromLength)) & !("Y" %in% names(ChromLength)))
                           {
                           Chromosomes = sort(as.numeric(names(ChromLength)))
                           }
                      ### Init vector with positions for chromosome endings on x-axis
                      Len = max(OverlayMat[OverlayMat$Chromosome==1,"Stop"])         
                      ### Init vector with positions for x-axis labeling                                 
                      LabelPos = Len/2
                      Lentmp = Len
                      chromosomes = Chromosomes[Chromosomes!=1]         
                      for(chromosome in chromosomes)
                                 {
                                 Lentmp = Lentmp + max(OverlayMat[OverlayMat$Chromosome==chromosome,"Stop"])
                                 OverlayMat[OverlayMat$Chromosome==chromosome,"Start"] = OverlayMat[OverlayMat$Chromosome==chromosome,"Start"]+max(Len)
                                 OverlayMat[OverlayMat$Chromosome==chromosome,"Stop"] = OverlayMat[OverlayMat$Chromosome==chromosome,"Stop"]+max(Len)
                                 labelpos = (max(OverlayMat[OverlayMat$Chromosome==chromosome,"Stop"])/2) + (max(Len)/2)
                                 LabelPos = c(LabelPos,labelpos)
                                 Len = c(Len,Lentmp)
                                 }
                      Len <<- Len
                      LabelPos <<- LabelPos           
                      OverlayMat <<- OverlayMat
                      Chromosomes <<- Chromosomes
                      }
                      

PlotOverlay = function(format="",Suffix="")
               {
               YcordsGain = OverlayMat$Gain
               YcordsDel = OverlayMat$Del
               Xcords = OverlayMat$Start
               Xcords = c(0,Xcords,0)
               YcordsGain = c(YcordsGain[1],YcordsGain,rev(YcordsGain)[1])
               YcordsDel = c(YcordsDel[1],YcordsDel,rev(YcordsDel)[1])    
               if(SummaryStat=="Mean")
                   {
                   Ylim=c(-6,4)
                   TickLabelPosition = seq(min(Ylim),max(Ylim),1)
                   TickLabel = seq(min(Ylim),max(Ylim),1)
                   Ylabelposition = 0
                   Ylabeltext = "Mean log2 copynumber change"
                   }
               if(SummaryStat=="Proportion")
                   {
                   Ylim=c(-1,1)
                   #Ylabel = "Proportion of samples"
                   TickLabelPosition = seq(-1,1,0.2)
                   TickLabel = abs(seq(-100,100,20))
                   Ylabelposition = c(-0.5,0.5) 
                   Ylabeltext = c("Proportion of samples (%) \n with deletions","Proportion of samples (%) \n with amplifications")
                   }
               if(format=="pdf")
                   {
                   pdf("CopyNumberOverlay.pdf",width=18,height=12,paper="special")
                   }
               par(las=1)            
               plot(1,xlim=c(0,max(OverlayMat$Stop)),ylim=Ylim,bty="n",xaxt="n",,yaxt="n",xlab="",ylab="",type="n")
               polygon(Xcords,YcordsGain,col="#B22222",border="#B22222")
               polygon(Xcords,YcordsDel,col="#4682B4",border="#4682B4")
               segments(Len,min(Ylim),Len,max(Ylim),lty=3,col="grey40",lwd=0.9)
               mtext(side=1,line=1,at=LabelPos,Chromosomes)
               mtext(side=2,line=1,las=3,at=Ylabelposition,text=Ylabeltext)
               axis(side=2,line=-3,at=TickLabelPosition,labels=TickLabel)
               if(format=="pdf")
                  {
                  dev.off()
                  }
               }



###########################
##### Main function  ######
###########################


RunOverlayAnalysis = function(Samples,species="Mouse",Method="Copywriter",resolution=20000,AberrationCutoff=0.2,SummaryStat="Proportion",ChromomsomesToRemove=c("X","Y"),format="pdf",Suffix="")
                         {
                         BinGenome(species=species,Method=Method,resolution=resolution,ChromomsomesToRemove=ChromomsomesToRemove)
                         SetOverlayMat(Samples)
                         fillOverlayMat(Samples,AberrationCutoff=AberrationCutoff,SummaryStat=SummaryStat)
                         print(ChromLength)
                         ReformatOverlayMat()
                         PlotOverlay()
                         }
                         
                         
                         
Samples = list.dirs(full.name=F,recursive=F)
Samples=Samples[-grep("CohortDatabase",Samples)]
Samples=Samples[-grep("misc",Samples)]
Samples=Samples[grep("PPT",Samples)]
pdf("Overlay.pdf")
RunOverlayAnalysis(Samples=Samples,SummaryStat="Proportion",ChromomsomesToRemove=c("X","Y"),format="pdf")
dev.off()
RunOverlayAnalysis(Samples=Samples,Method="HMMCopy",resolution=20000)