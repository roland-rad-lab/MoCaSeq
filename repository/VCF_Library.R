#!/usr/bin/env Rscript



bedpe_from_jabba <- function (data)
{

	data_aux <- data %>%
		dplyr::mutate (id_sample=paste (id,sample,sep="__"),info.mateid_sample=paste (info.mateid,sample,sep="__"),info.event_sample=paste (info.event,sample,sep="__")) %>%
		data.frame

	inversion_head_to_head <- data_aux %>%
		dplyr::filter (info.svtype=="BND",stringr::str_detect (alt,"^\\]|^\\[")) %>%
		dplyr::inner_join (data_aux,by=c("id_sample"="info.mateid_sample"),suffix=c("1","2")) %>%
		dplyr::mutate (sample=sample1,SVtype="h2hINV",strand1="+",strand2="+") %>%
		dplyr::mutate (start1=pos1-1,start2=pos2-1,name=paste (id1,id2,sep=""),score=1) %>%
		dplyr::rename (end1=pos1,end2=pos2) %>%
		dplyr::select (sample,chrom1,start1,end1,chrom2,start2,end2,name,score,strand1,strand2) %>%
		data.frame

	inversion_tail_to_tail <- data_aux %>%
		dplyr::filter (info.svtype=="BND",stringr::str_detect (alt,"\\]$|\\[$")) %>%
		dplyr::inner_join (data_aux,by=c("id_sample"="info.mateid_sample"),suffix=c("1","2")) %>%
		dplyr::mutate (sample=sample1,SVtype="t2tINV",strand1="-",strand2="-") %>%
		dplyr::mutate (start1=pos1-1,start2=pos2-1,name=paste (id1,id2,sep=""),score=1) %>%
		dplyr::rename (end1=pos1,end2=pos2) %>%
		dplyr::select (sample,chrom1,start1,end1,chrom2,start2,end2,name,score,strand1,strand2) %>%
		data.frame

	result <- rbind (inversion_head_to_head,inversion_tail_to_tail)

	return (result)
}

bedpe_from_manta <- function (data)
{
	data_aux <- data %>%
		dplyr::filter (alt!="<INS>") %>%
		dplyr::mutate (info.end=as.numeric (info.end)) %>%
		dplyr::mutate (id_sample=paste (id,sample,sep="__"),info.mateid_sample=paste (info.mateid,sample,sep="__"),info.event_sample=paste (info.event,sample,sep="__")) %>%
		data.frame

	deletions <- data_aux %>%
		dplyr::filter (info.svtype=="DEL") %>%
		dplyr::mutate (chrom1=chrom,end1=pos,chrom2=chrom,end2=info.end,name=paste (id,info.svtype,sep="__"),score=1,strand1="+",strand2="-") %>%
		dplyr::mutate (start1=pos-1,start2=info.end-1) %>%
		dplyr::select (sample,chrom1,start1,end1,chrom2,start2,end2,name,score,strand1,strand2) %>%
		data.frame

	duplications <- data_aux %>%
		dplyr::filter (info.svtype=="DUP") %>%
		dplyr::mutate (chrom1=chrom,end1=pos,chrom2=chrom,end2=info.end,name=paste (id,info.svtype,sep="__"),score=1,strand1="-",strand2="+") %>%
		dplyr::mutate (start1=pos-1,start2=info.end-1) %>%
		dplyr::select (sample,chrom1,start1,end1,chrom2,start2,end2,name,score,strand1,strand2) %>%
		data.frame

	# Inversions are described in terms of breakends (type BND) according to the VCF 4.2 spec

	inversion_head_to_head <- data_aux %>%
		dplyr::filter (info.svtype=="BND",stringr::str_detect (alt,"^\\]|^\\[")) %>%
		dplyr::inner_join (data_aux,by=c("id_sample"="info.mateid_sample"),suffix=c("1","2")) %>%
		dplyr::mutate (sample=sample1,SVtype="h2hINV",strand1="+",strand2="+") %>%
		dplyr::mutate (start1=pos1-1,end1=pos1,start2=pos2-1,end2=pos2,name=paste (id1,id2,SVtype,sep="__"),score=1) %>%
		dplyr::select (sample,chrom1,start1,end1,chrom2,start2,end2,name,score,strand1,strand2) %>%
		data.frame

	inversion_tail_to_tail <- data_aux %>%
		dplyr::filter (info.svtype=="BND",stringr::str_detect (alt,"\\]$|\\[$")) %>%
		dplyr::inner_join (data_aux,by=c("id_sample"="info.mateid_sample"),suffix=c("1","2")) %>%
		dplyr::mutate (sample=sample1,SVtype="t2tINV",strand1="-",strand2="-") %>%
		dplyr::mutate (start1=pos1-1,end1=pos1,start2=pos2-1,end2=pos2,name=paste (id1,id2,SVtype,sep="__"),score=1) %>%
		dplyr::select (sample,chrom1,start1,end1,chrom2,start2,end2,name,score,strand1,strand2) %>%
		data.frame

	result <- rbind (deletions,duplications,inversion_head_to_head,inversion_tail_to_tail)

	return (result)
}

