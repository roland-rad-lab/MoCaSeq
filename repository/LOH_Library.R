#!/usr/bin/Rscript

##########################################################################################
##
## LOH_Library.R
##
## Various functions needed for LOH_GenerateVariantTable.R
##
##########################################################################################

library("GenomicRanges")
library("Biostrings")

LOH_FilterReads = function(data)
{
	data=data[data[,"Chrom"] %in% c(1:22,"X","Y"),]
	data=data[-grep(",",data[,"Alt"]),]
	data[,"MapQ"]=as.numeric(as.character(data[,"MapQ"]))
	data[,"Frequency"]=as.numeric(as.character(data[,"Frequency"]))
	data=data[data[,"RefCount"] + data[,"AltCount"] >= 10,]
	data=data[data[,"MapQ"]>=60,]
	return(data)
}

LOH_MergeVariants = function(tumor,normal)
{
	tumor[,"UniquePos"]=paste(tumor[,"Chrom"], tumor[,"Pos"], tumor[,"Ref"], tumor[,"Alt"], sep="_")
	normal[,"UniquePos"]=paste(normal[,"Chrom"], normal[,"Pos"], normal[,"Ref"], normal[,"Alt"], sep="_")
	variants <- merge(x=normal,y=tumor,by="UniquePos",all=F)
	variants = variants[,c("UniquePos", "Chrom.x", "Pos.x", "Ref.x", "Alt.x", "Frequency.x", "RefCount.x", "AltCount.x", "Frequency.y","RefCount.y", "AltCount.y")]
	colnames(variants) = c("UniquePos", "Chrom", "Pos", "Ref", "Alt", "Normal_Freq", "Normal_RefCount", "Normal_AltCount", "Tumor_Freq", "Tumor_RefCount", "Tumor_AltCount")
	variants[,"Chrom"]=as.character(variants[,"Chrom"])
	return(variants)
}

LOH_DefineDictionaries = function()
{
	dicts=list()
	dict_for=c("AC","AG","CA","GA")
	dict_rev=c("TC","TG","CT","GT")
	dict_unamb=c(dict_for,dict_rev)

	dicts[["dict_for"]]=dict_for
	dicts[["dict_rev"]]=dict_rev
	dicts[["dict_unamb"]]=dict_unamb
	return(dicts)
}

LOH_AssignStatus = function(variants,dict_unamb)
{
	results=list()

	variants[,"Genotype"]=paste(variants[,"Ref"], variants[,"Alt"],sep="")
	amb=variants[!variants[,"Genotype"] %in% dict_unamb,]
	indel=amb[!amb[,"Ref"] %in% c("T","C","A","G") & amb[,"Alt"] %in% c("T","C","A","G"),]
	amb=amb[amb[,"Ref"] %in% c("T","C","A","G") & amb[,"Alt"] %in% c("T","C","A","G"),]
	unamb=variants[variants[,"Genotype"] %in% dict_unamb,]

	unamb[,"Plot_Freq"] = 0
	unamb[unamb[,"Genotype"] %in% dict_for,"Plot_Freq"]=unamb[unamb[,"Genotype"] %in% dict_for,"Tumor_Freq"]
	unamb[unamb[,"Genotype"] %in% dict_rev,"Plot_Freq"]=1-unamb[unamb[,"Genotype"] %in% dict_rev,"Tumor_Freq"]
	unamb["strand_vector"]="NA"

	#indels are not assigned to A- or B-Allele but always assigned to the B-Allele
	indel[,"Plot_Freq"] = indel[,"Tumor_Freq"]
	indel["strand_vector"]="NA"

	results[["amb"]]=amb
	results[["indel"]]=indel
	results[["unamb"]]=unamb
	return(results)
}

LOH_GenerateLookupTable = function(name,amb)
{
	export=data.frame(amb[,c("Chrom","Pos")])
	export[,"PosStart"]=export[,"Pos"]-100
	export[,"PosEnd"]=export[,"Pos"]+100
	export$export=paste(export[,"Chrom"],":",export[,"PosStart"], "-",export[,"PosEnd"],sep="")

	write.table(export$export,paste(name,".lookup.tab",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
}

LOH_ExtractSurroundingNucleotides = function(name,genome_fasta)
{
	system(paste("parallel -a ",name,".lookup.tab 'samtools faidx ",genome_fasta," {}' > ",name,".found.tab",sep=""))
}

LOH_ImportSurroundingNucleotides = function(name)
{
	import=data.frame(Seq=unstrsplit(paste(name,".found.tab",sep="")))
	import[,"Pos"]=row.names(import)
	row.names(import)=NULL
	import=import[,c(2,1)]

	system(paste("rm ",name,".lookup.tab",sep=""))
	system(paste("rm ",name,".found.tab",sep=""))

	return(import)
}

LOH_AssignAlleles = function(amb,import,dict_unamb)
{
	strand_vector=c()

	for (i in 1:nrow(import))
	{
		j=1
		while (!paste(substr(import[i,2],101-j,101-j),substr(import[i,2],101+j,101+j),sep="") %in% dict_unamb & j<100)
		{
			j=j+1
		}
		if (paste(substr(import[i,2],11-j,11-j),substr(import[i,2],11+j,11+j),sep="") %in% c("AC","AG", "TC","TG"))
		{
			strand_vector=c(strand_vector,"TOP")
		} else {
			strand_vector=c(strand_vector,"BOTTOM")
		}
	}

	amb$strand_vector=strand_vector

	for (i in 1:nrow(amb))
	{
		if (amb[i,"strand_vector"]=="TOP" & amb[i,"Ref"] %in% c("A","C"))
		{
			amb[i,"Plot_Freq"]=1-amb[i,"Tumor_Freq"]
		}
		if (amb[i,"strand_vector"]=="TOP" & amb[i,"Ref"] %in% c("T","G"))
		{
			amb[i,"Plot_Freq"]=amb[i,"Tumor_Freq"]
		}
		if (amb[i,"strand_vector"]=="BOTTOM" & amb[i,"Ref"] %in% c("A","C"))
		{
			amb[i,"Plot_Freq"]=amb[i,"Tumor_Freq"]
		}
		if (amb[i,"strand_vector"]=="BOTTOM" & amb[i,"Ref"] %in% c("T","G"))
		{
			amb[i,"Plot_Freq"]=1-amb[i,"Tumor_Freq"]
		}
	}

	return(amb)
}