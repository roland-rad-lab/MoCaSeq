#!/usr/bin/Rscript

##########################################################################################
##
## Chromothripsis_SimulateCopyNumberStates.R
##
## Simulate a progressive tumor model in a diploid organism given a list of structural rearrangements.
## The chromosomes are simulated as lists of chromosomal ranges which are modified by the rearrangements.
##
##########################################################################################

message("\n###Simulate Progressive Tumor###")
options(warn=-1)
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-i", "--input"),type="character",default=NULL,help="rearrangement list"),
  make_option(c("-o", "--organism"),type="character",default=NULL,help="human|mouse"),
  make_option(c("-c", "--chrom"),type="character",default=NULL,help="chromosome to be modeled"),
  make_option(c("-s", "--steps"),type="integer",default=30,help="number of MC simulations [default = %default]"),
  make_option(c("-a", "--abort"),type="integer",default=1000,help="number of tries, before simulation is aborted [default = %default]"),
  make_option(c("-v", "--verbose"),type="integer",default=0,help="verbose level (0-3) [default = %default]"),
  make_option(c("-n", "--name"),type="character",default="Sample1",help="sample name [default = %default]"),
  make_option(c("-f", "--format"),type="character",default="tif",help="output format (tif|emf) [default = %default]")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

human_chroms <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,
                  135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,
                  50818468,156040895,57227415)
names(human_chroms)=c(1:22,"X","Y")
mouse_chroms <- c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,
                  122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566,171031299,91744698)
names(mouse_chroms)=c(1:19,"X","Y")


if(is.null(opt$input)){
  print_help(opt_parser)
  stop("You have to specify a rearrangement list.", call.=FALSE)
} else{
  tryCatch({
    rearrangementList <- read.csv(opt$input,header=T,sep="\t",stringsAsFactors=F)
  },
  error = function(e){
    stop("Error reading rearrangement list.")
  })
}
if(is.null(opt$organism)){
  print_help(opt_parser)
  stop("You have to specify an organism (either human or mouse).", call.=FALSE)
} else{
  if(is.null(opt$chrom)){
    print_help(opt_parser)
    stop("You have to specify a chromosome to model.", call.=FALSE)
  } else {
    if(opt$organism == "human" & opt$chrom %in% names(human_chroms)) chromlength <- human_chroms[opt$chrom]
    else if(opt$organism == "mouse" & opt$chrom %in% names(mouse_chroms)) chromlength <- mouse_chroms[opt$chrom]
    else{
      print_help(opt_parser)
      stop("The given organism and the specified chromsome do not match.", call.=FALSE)
    }
    rearrangementList = rearrangementList[rearrangementList$donChr==opt$chrom & rearrangementList$accChr==opt$chrom,]
  }
}
if(nrow(rearrangementList)<4) stop("There are too few (< 4) rearrangements on the specified chromosome. No simulation possible.",call.=FALSE)
if(!opt$format %in% c("tif","emf")) stop("The output format has to be either tif or emf.",call.=FALSE)
MC_steps <- opt$steps # default = 100
verboseLevel = opt$verbose # this variable controls how much information is displayed on the command line; default = 0


suppressPackageStartupMessages(library(devEMF))
suppressPackageStartupMessages(library(GenomicRanges))
start_chrom <- data.frame(start=1,end=chromlength)

# auxiliary functions
# auxiliary function to check whether the insertion sites still exist in at least the simulated chromosomes;
# if this is the case, return the chromosomal positions, otherwise return NULL
findIndices <- function(start,end,chroms){
	chr_available <- c()

	if(computeCN(start,chroms[[1]]) > 0 & computeCN(end,chroms[[1]]) > 0) chr_available <- c(chr_available,1)
	if(computeCN(start,chroms[[2]]) > 0 & computeCN(end,chroms[[2]]) > 0) chr_available <- c(chr_available,2)

	if(is.null(chr_available)) return(NULL)

	currentChr <- chr_available[sample(1:length(chr_available),1)]

	startInd <- which((chroms[[currentChr]][,1] <= start & chroms[[currentChr]][,2] >= start) | (chroms[[currentChr]][,2] <= start & chroms[[currentChr]][,1] >= start))
	startInd <- startInd[sample(1:length(startInd),1)]
	endInd <- which((chroms[[currentChr]][,1] <= end & chroms[[currentChr]][,2] >= end) | (chroms[[currentChr]][,2] <= end & chroms[[currentChr]][,1] >= end))
	endInd <- endInd[sample(1:length(startInd),1)]

	if(startInd < endInd){
		return(c(startInd,endInd,currentChr,start,end))
	}
	else if(startInd == endInd){
		if(chroms[[currentChr]][startInd,2] < chroms[[currentChr]][startInd,1]) return(c(startInd,endInd,currentChr,end,start))
		else return(c(startInd,endInd,currentChr,start,end))
	}
	else{
		return(c(endInd,startInd,currentChr,end,start))
	}
}

# auxiliary function for insertion of deletions
insertDel <- function(position,start,end,chrom){
	# reassemble chromosome
	newChrom <- chrom[1:position[1]-1,]
	newChrom <- rbind(newChrom,data.frame(start=chrom[position[1],1],end=start))
	newChrom <- rbind(newChrom,data.frame(start=end,end=chrom[position[2],2]))
	newChrom <- rbind(newChrom,chrom[position[2]+1:nrow(chrom),])
	newChrom <- newChrom[!is.na(newChrom[,1]),]

	return(newChrom)
}

# auxiliary function for insertion of duplications
insertDup <- function(position,start,end,chrom){
	# reassemble chromosome
	newChrom <- chrom[1:position[2]-1,]
	newChrom <- rbind(newChrom,data.frame(start=chrom[position[2],1],end=end))
	newChrom <- rbind(newChrom,data.frame(start=start,end=chrom[position[1],2]))
	newChrom <- rbind(newChrom,chrom[position[1]+1:nrow(chrom),])
	newChrom <- newChrom[!is.na(newChrom[,1]),]

	return(newChrom)
}

# auxiliary function for insertion of inversions
insertInv <- function(position,start,end,chrom){
	startRev = 1
	if(chrom[position[1],2] < chrom[position[1],1]) startRev = -1
	endRev = 1
	if(chrom[position[2],2] < chrom[position[2],1]) endRev = -1

	# reassemble chromosome
	newChrom <- chrom[1:position[1]-1,]
	newChrom <- rbind(newChrom,data.frame(start=chrom[position[1],1],end=start))
	if(position[1] == position[2]){
		newChrom <- rbind(newChrom,data.frame(start=end,end=start+startRev))
	}
	else{
		newChrom <- rbind(newChrom,data.frame(start=end,end=chrom[position[2],1]))
		if(position[2] - position[1] > 1){
			for(i in 1:(position[2]-position[1]-1)){
				newChrom <- rbind(newChrom,data.frame(start=chrom[position[2]-i,2],end=chrom[position[2]-i,1]))
			}
		}
		if(chrom[position[1],2] != start){
			newChrom <- rbind(newChrom,data.frame(start=chrom[position[1],2],end=start+startRev))
		}
	}
	if(chrom[position[2],2] != end){
		newChrom <- rbind(newChrom,data.frame(start=end+endRev,end=chrom[position[2],2]))
	}
	newChrom <- rbind(newChrom,chrom[position[2]+1:nrow(chrom),])
	newChrom <- newChrom[!is.na(newChrom[,1]),]

	return(newChrom)
}

# auxiliary function which tries to combine neighbouring fragments from simulated chromosomes
currateChromosome <- function(chrom){
	newChrom <- chrom[1,]

	for(i in 2:nrow(chrom)){
		if((chrom[i-1,1] <= chrom[i-1,2] & chrom[i,1] <= chrom[i,2] & chrom[i,1] == chrom[i-1,2] + 1) | (chrom[i-1,2] <= chrom[i-1,1] & chrom[i,2] <= chrom[i,1] & chrom[i,1] == chrom[i-1,2] - 1)){
			newChrom <- rbind(newChrom[0:(nrow(newChrom)-1),],data.frame(start=chrom[i-1,1],end=chrom[i,2]))
		}
		else{
			newChrom <- rbind(newChrom,chrom[i,])
		}
	}

	return(newChrom)
}

# compute the copynumber state in the simulated chromosome pair at a given position
computeCN <- function(pos,chrom){
	currentState = 0

	for(i in 1:nrow(chrom)){
		if((chrom[i,1] <= pos & chrom[i,2] >= pos) | (chrom[i,2] <= pos & chrom[i,1] >= pos)) currentState = currentState + 1
	}

	return(currentState)
}

# compute the number of different copynumber states for a simulated pair of chromosomes
computeMaxCN <- function(chroms){
	outlist <- rbind(chroms[[1]],chroms[[2]])
	for(i in 1:nrow(outlist)){
		if(outlist[i,2] < outlist[i,1]){
			a = outlist[i,1]
			outlist[i,1] = outlist[i,2]
			outlist[i,2] = a
		}
	}
	outlist$chr = paste0("chr",opt$chrom)
	outlist = makeGRangesFromDataFrame(outlist)

	return(length(unique(runValue(coverage(outlist)[[paste0("chr",opt$chrom)]]))))
}

# visualize simulated chromosomes as copy number plots
visualizeChrom <- function(chroms,file){
	outlist <- rbind(chroms[[1]],chroms[[2]])

	for(i in 1:nrow(outlist)){
		if(outlist[i,2] < outlist[i,1]){
			a = outlist[i,1]
			outlist[i,1] = outlist[i,2]
			outlist[i,2] = a
		}
	}
	outlist$chr = paste0("chr",opt$chrom)
	outlist = makeGRangesFromDataFrame(outlist)

	cov_pos = runLength(coverage(outlist)[[paste0("chr",opt$chrom)]])
	cov_val = runValue(coverage(outlist)[[paste0("chr",opt$chrom)]])
	CNstate = data.frame(start=1,end=cumsum(cov_pos),value=cov_val)
	for(i in 2:nrow(CNstate)) CNstate$start[i] = CNstate$end[i-1]+1
	CNstate$value[CNstate$value==0] = 0.5
	CNstate$value = log2(CNstate$value/2)

	tiff(file,1855,456)
		plot(x=0,y=0,xlim=c(0,ceiling(chromlength/10000000)*10000000),ylim=c(round(min(CNstate$value))-1,round(max(CNstate$value))+1),
		     xlab='Genomic Position',ylab='log2 CN State',type='n',frame.plot=F)
		for(i in 1:nrow(CNstate)){
			xpos <- seq(CNstate$start[i],CNstate$end[i],by=1000)
			points(xpos,rnorm(length(xpos),mean=CNstate$value[i],sd=0.1),pch='.')
			segments(CNstate$start[i],CNstate$value[i],CNstate$end[i],CNstate$value[i],lwd=2,col='red')
		}
	garbage <- dev.off()
}


# main function
#
# in each step, the algorithm tries to insert i randomly chosen rearrangements from the input list, where i goes from one
# to the total number of rearragements; this procedure is repeated <MC_steps> times for every i; if a rearrangement can be
# inserted at more than one site, one of them is chosen randomly; if a rearrangement cannot be inserted (because at least
# one of the corresponding positions doesn't exist anymore), it is discarded and a new one is chosen randomly; if there are
# no rearrangements left in the input list and it was not possible to insert i rearrangements, the procedure is repeated
# until <MC_steps> rearrangement cycles have been performed successfully

CNstatesAll <- setNames(data.frame(matrix(ncol=MC_steps,nrow=0)),paste0("step_",c(1:MC_steps)))
chromosomes <- list() # a list storing all chromosomal rearrangements (for debugging)
savePath = paste0(opt$name,"/results/Chromothripsis/Chr",opt$chrom,"/",opt$name,"_simulation_chr",opt$chrom)
system(paste0("mkdir -p ",savePath))
stop_flag = F

if(verboseLevel>0) message(paste0(nrow(rearrangementList)," rearrangements to simulate.\n"))

for(i in 1:nrow(rearrangementList)){
	# i stores the number of rearrangemets to insert in the current step
  if(verboseLevel>0) message(paste0("\nSimulating ",i," rearrangement(s)."))

	CNstates <- c() # this vector holds the numbers of different copy number states for the current rearrangement cycle
	j = 1 # j counts the number of successful rearrangement cycles of the current loop
	f = 0 # f counts the number of failed rearrangement cycles of the current loop

	while(j <= MC_steps){
  	if(verboseLevel>1) message(paste0('\tDoing ',j+f,'. simulation with ',i,' rearrangement(s).'))

	  # permute the rearrangement list to produce the random insertion ordera and initialize chromosomes
		currentEvents <- rearrangementList[sample(1:nrow(rearrangementList)),]
		chroms <- list(chrom1=start_chrom,chrom2=start_chrom)

		k = 1 # counts the number of successfull rearrangements in the current cycle
		e = 0 # counts the number of failed rearrangements in the current cycle

		while(k <= i & k + e <= nrow(currentEvents)){
			# find the start and end position of the current rearrangement and check whether these positions still exist in the simulated chromosomes
		  start <- currentEvents$donPos[k+e]
			end <- currentEvents$accPos[k+e]
			pos <- findIndices(start,end,chroms)

			if(!is.null(pos)){ # do the following only if the insertion sites exist
				if(verboseLevel>2) message(paste0('\t\tInserting ',k,'. rearrangement (',k+e,'. in list).'))

			  # perform the actual rearrangement
			  switch(as.character(currentEvents$Type[k+e]),
					'3to5' = {chroms[[pos[3]]] <- insertDel(pos[1:2],pos[4],pos[5],chroms[[pos[3]]])},
					'5to3' = {chroms[[pos[3]]] <- insertDup(pos[1:2],pos[4],pos[5],chroms[[pos[3]]])},
					{chroms[[pos[3]]] <- insertInv(pos[1:2],pos[4],pos[5],chroms[[pos[3]]])}
				)
				chroms[[pos[3]]] <- currateChromosome(chroms[[pos[3]]]) # try to compact the chromosomal ranges
				k = k + 1 # increase the counter for successfull insertions
			}
			else{ # if it is not possible to insert the current rearrangement, increase the error counter
				if(verboseLevel>2) message(paste0('\t\tCouldn\'t insert ',k,'. rearrangement (',k+e,'. in list). Trying next.'))
				e = e + 1
			}
		}

		# if the desired number of rearrangements (i) could be successfully inserted, store the number of different copynumber
	  # states and the chromosome range list and increase the counter for successful rearrangement cycles...
		if(k == i + 1){
		  if(verboseLevel==1 & j%%10==0) message(paste0('\t',j,' simulation(s) done.'))
		  if(verboseLevel>1) message(paste0('\t',j,'. success.'))
		  chromosomes[[paste0(i,'_',j)]] = chroms
			CNstates <- c(CNstates,computeMaxCN(chroms))
			j = j + 1
		}
		else{ # ...otherwise increase the rearrangement cycle error count or abort
		  if(f+1 >= opt$abort){ # this break condition is used to stop the simulation, if it is not possible to insert all rearrangements after a given number of tries
		    stop_flag = T
		    message(paste0('\t',f+1,'. fail. Aborting simulation.'))
		    break
		  } else{
		    if(verboseLevel>1) message(paste0('\t',f+1,'. fail. Trying again.'))
		    f = f + 1
		  }
		}
	}
	if(stop_flag) break

	CNstatesAll <- rbind(CNstatesAll,CNstates)
	colnames(CNstatesAll) = paste0("step_",c(1:MC_steps))
	save(CNstatesAll,file=paste0(savePath,"/CNstatesAll.RData"))
	save(chromosomes,file=paste0(savePath,"/chromosomes.RData"))

	# generate copynumber plot for the chromosome with the highest CN
	visualizeChrom(chromosomes[[paste0(i,"_",which(CNstatesAll[i,]==max(CNstatesAll[i,]))[1])]],paste0(savePath,"/Profile_insert",i,".tif"))
}

# Plotting routine
if(opt$format=="tif")
   {
    tiff(paste0(opt$name,"/results/Chromothripsis/Chr",opt$chrom,"/",opt$name,".chr",opt$chrom,".CopyNumberSimulation.tif"),1600,1600,res=200)
    }
if(opt$format=="emf")
     {
     emf(paste0(opt$name,"/results/Chromothripsis/Chr",opt$chrom,"/",opt$name,".chr",opt$chrom,".CopyNumberSimulation.emf"),bg="white", width=8,height=8,coordDPI = 200)
     }
par(oma=c(4,4,2,2))
x_size = nrow(rearrangementList)+(-nrow(rearrangementList)%%10)
y_size = max(apply(CNstatesAll,1,mean)+1.96*apply(CNstatesAll,1,sd)/sqrt(MC_steps))+(-max(apply(CNstatesAll,1,mean)+1.96*apply(CNstatesAll,1,sd)/sqrt(MC_steps))%%5)
plot(x=0,y=0,las=1,frame.plot=F,xaxt='n',yaxt='n',xlim=c(0,x_size),ylim=c(-2,y_size),xlab='',ylab='',type='n')

points(apply(CNstatesAll,1,mean),pch=20,cex=1.2)
axis(side=2,las=1,at=seq(0,y_size,by=5),labels=seq(0,y_size,by=5),lwd=2.2,cex.axis=2)
mtext(side=2,las=3,at=y_size/2,line=3,c('Number of copy number states'),cex=2)
axis(side=1,at=seq(0,x_size,by=10),labels=seq(0,x_size,by=10),line=-4.5,lwd=2.2,cex.axis=2)
mtext(side=1,las=1,at=x_size/2,line=-1,c('Number of breakpoints'),cex=2)

segments(c(1:nrow(CNstatesAll)),apply(CNstatesAll,1,mean)-1.96*apply(CNstatesAll,1,sd)/sqrt(MC_steps),c(1:nrow(CNstatesAll)),apply(CNstatesAll,1,mean)+1.96*apply(CNstatesAll,1,sd)/sqrt(MC_steps),cex=2)
#points(x=nrow(rearrangementList),y=sample_pos,pch=23,col='darkred',bg='darkred')
#text(x=nrow(rearrangementList),y=sample_pos,labels=opt$name,pos=1)
text(3,y_size,paste("n=",nrow(rearrangementList),sep=""),cex=2)

print("Done")
garbage <- dev.off()