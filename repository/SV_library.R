#!/usr/bin/env Rscript


sv_from_manta <- function (data)
{
	manta_data_aux <- data %>%
		filter (alt!="<INS>") %>%
		mutate (info.end=as.numeric (info.end)) %>%
		mutate (id_sample=paste (id,sample,sep="__"),info.mateid_sample=paste (info.mateid,sample,sep="__"),info.event_sample=paste (info.event,sample,sep="__")) %>%
		data.frame

	#print (head (manta_data_aux))

	deletions <- manta_data_aux %>%
		filter (info.svtype=="DEL") %>%
		mutate (chrom1=chrom,pos1=pos,chrom2=chrom,pos2=info.end,SVtype=info.svtype,strand1="+",strand2="-") %>%
		dplyr::select (sample,chrom1,pos1,chrom2,pos2,SVtype,strand1,strand2) %>%
		data.frame

	duplications <- manta_data_aux %>%
		filter (info.svtype=="DUP") %>%
		mutate (chrom1=chrom,pos1=pos,chrom2=chrom,pos2=info.end,SVtype=info.svtype,strand1="-",strand2="+") %>%
		dplyr::select (sample,chrom1,pos1,chrom2,pos2,SVtype,strand1,strand2) %>%
		data.frame

	# Inversions are described in terms of breakends (type BND) according to the VCF 4.2 spec

	#print (manta_data_aux %>% filter (sample=="Normal",info.event=="MantaBND:84305:2:4:0:0:0:0") %>% data.frame)
	#print (manta_data_aux %>% filter (sample=="Normal",info.event=="MantaBND:450650:0:1:0:0:0:0") %>% data.frame)

	inversion_head_to_head <- manta_data_aux %>%
		filter (info.svtype=="BND",stringr::str_detect (alt,"^\\]|^\\[")) %>%
		inner_join (manta_data_aux,by=c("id_sample"="info.mateid_sample"),suffix=c("1","2")) %>%
		mutate (sample=sample1,SVtype="h2hINV",strand1="+",strand2="+") %>%
		dplyr::select (sample,chrom1,pos1,chrom2,pos2,SVtype,strand1,strand2) %>%
		data.frame

	inversion_tail_to_tail <- manta_data_aux %>%
		filter (info.svtype=="BND",stringr::str_detect (alt,"\\]$|\\[$")) %>%
		inner_join (manta_data_aux,by=c("id_sample"="info.mateid_sample"),suffix=c("1","2")) %>%
		mutate (sample=sample1,SVtype="t2tINV",strand1="-",strand2="-") %>%
		dplyr::select (sample,chrom1,pos1,chrom2,pos2,SVtype,strand1,strand2) %>%
		data.frame

	result <- rbind (deletions,duplications,inversion_head_to_head,inversion_tail_to_tail)

	return (result)
}

sv_from_delly <- function (data)
{
	delly_data_aux <- data %>%
		filter (alt!="<INS>") %>%
		mutate (info.end=as.numeric (info.end)) %>%
		mutate (id_sample=paste (id,sample,sep="__")) %>%
		data.frame

	deletions <- delly_data_aux %>%
		filter (info.svtype=="DEL") %>%
		mutate (chrom1=chrom,pos1=pos,chrom2=chrom,pos2=info.end,SVtype=info.svtype,strand1="+",strand2="-") %>%
		dplyr::select (sample,chrom1,pos1,chrom2,pos2,SVtype,strand1,strand2) %>%
		data.frame

	duplications <- delly_data_aux %>%
		filter (info.svtype=="DUP") %>%
		mutate (chrom1=chrom,pos1=pos,chrom2=chrom,pos2=info.end,SVtype=info.svtype,strand1="-",strand2="+") %>%
		dplyr::select (sample,chrom1,pos1,chrom2,pos2,SVtype,strand1,strand2) %>%
		data.frame

	inversion_head_to_head <- delly_data_aux %>%
		filter (info.svtype=="BND",stringr::str_detect (alt,"^\\]|^\\[")) %>%
		mutate (chrom1=chrom,pos1=pos,chrom2=info.chr2,pos2=info.pos2,SVtype="h2hINV",strand1="+",strand2="+") %>%
		dplyr::select (sample,chrom1,pos1,chrom2,pos2,SVtype,strand1,strand2) %>%
		data.frame

	inversion_tail_to_tail <- delly_data_aux %>%
		filter (info.svtype=="BND",stringr::str_detect (alt,"\\]$|\\[$")) %>%
		mutate (chrom1=chrom,pos1=pos,chrom2=info.chr2,pos2=info.pos2,SVtype="t2tINV",strand1="-",strand2="-") %>%
		dplyr::select (sample,chrom1,pos1,chrom2,pos2,SVtype,strand1,strand2) %>%
		data.frame

	result <- rbind (deletions,duplications,inversion_head_to_head,inversion_tail_to_tail)

	return (result)
}


