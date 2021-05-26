#!/usr/bin/env Rscript

# The manta VCF is somatic variants only
sv_from_manta <- function (data)
{
	manta_data_aux <- data %>%
		bind_cols (setNames(as.data.frame(apply(stringr::str_split_fixed (.$pr,",",2),2,as.numeric),stringsAsFactors=F),c("pr.ref","pr.alt"))) %>%
		bind_cols (setNames(as.data.frame(apply(stringr::str_split_fixed (.$sr,",",2),2,as.numeric),stringsAsFactors=F),c("sr.ref","sr.alt"))) %>%
		filter (alt!="<INS>") %>%
		mutate (info.end=as.numeric (info.end),info.svlen=as.numeric (info.svlen)) %>%
		mutate (id_sample=paste (id,sample,sep="__"),info.mateid_sample=paste (info.mateid,sample,sep="__"),info.event_sample=paste (info.event,sample,sep="__")) %>%
		data.frame

	#print (head (manta_data_aux))

	deletions <- manta_data_aux %>%
		filter (info.svtype=="DEL") %>%
		mutate (chrom1=chrom,pos1=pos,chrom2=chrom,pos2=info.end,SVtype=info.svtype,strand1="+",strand2="-") %>%
		mutate (gt=case_when (
				      pr.ref + sr.ref == 0 ~ "1/1",
				      pr.alt + sr.alt == 0 ~ "0/0",
				      TRUE ~ "0/1"
				      )) %>%
		dplyr::select (sample,chrom1,pos1,chrom2,pos2,SVtype,strand1,strand2,gt) %>%
		data.frame

	duplications <- manta_data_aux %>%
		filter (info.svtype=="DUP") %>%
		mutate (chrom1=chrom,pos1=pos,chrom2=chrom,pos2=info.end,SVtype=info.svtype,strand1="-",strand2="+") %>%
		mutate (gt=case_when (
				      pr.ref + sr.ref == 0 ~ "1/1",
				      pr.alt + sr.alt == 0 ~ "0/0",
				      TRUE ~ "0/1"
				      )) %>%
		dplyr::select (sample,chrom1,pos1,chrom2,pos2,SVtype,strand1,strand2,gt) %>%
		data.frame

	# Inversions are described in terms of breakends (type BND) according to the VCF 4.2 spec

	#print (manta_data_aux %>% filter (sample=="Normal",info.event=="MantaBND:84305:2:4:0:0:0:0") %>% data.frame)
	#print (manta_data_aux %>% filter (sample=="Normal",info.event=="MantaBND:450650:0:1:0:0:0:0") %>% data.frame)

	inversion_head_to_head <- manta_data_aux %>%
		filter (info.svtype=="BND",stringr::str_detect (alt,"^\\]|^\\[")) %>%
		inner_join (manta_data_aux,by=c("id_sample"="info.mateid_sample"),suffix=c("1","2")) %>%
		mutate (sample=sample1,SVtype="h2hINV",strand1="+",strand2="+") %>%
		mutate (gt=case_when (
				      pr.ref1 + pr.ref2 + sr.ref1 + sr.ref2 == 0 ~ "1/1",
				      pr.alt1 + pr.alt2 + sr.alt1 + sr.alt2 == 0 ~ "0/0",
				      TRUE ~ "0/1"
				      )) %>%
		dplyr::select (sample,chrom1,pos1,chrom2,pos2,SVtype,strand1,strand2,gt) %>%
		data.frame

	inversion_tail_to_tail <- manta_data_aux %>%
		filter (info.svtype=="BND",stringr::str_detect (alt,"\\]$|\\[$")) %>%
		inner_join (manta_data_aux,by=c("id_sample"="info.mateid_sample"),suffix=c("1","2")) %>%
		mutate (sample=sample1,SVtype="t2tINV",strand1="-",strand2="-") %>%
mutate (gt=case_when (
				      pr.ref1 + pr.ref2 + sr.ref1 + sr.ref2 == 0 ~ "1/1",
				      pr.alt1 + pr.alt2 + sr.alt1 + sr.alt2 == 0 ~ "0/0",
				      TRUE ~ "0/1"
				      )) %>%

		dplyr::select (sample,chrom1,pos1,chrom2,pos2,SVtype,strand1,strand2,gt) %>%
		data.frame

	result <- rbind (deletions,duplications,inversion_head_to_head,inversion_tail_to_tail)

	return (result)
}

sv_from_delly <- function (data)
{
	delly_data_aux <- data %>%
		#filter (ft=="PASS") %>%# The FILTER is PASS but FORMAT/FT is LowQual for all samples??
		filter (alt!="<INS>") %>%
		mutate (info.end=as.numeric (info.end),info.svlen=as.numeric (info.svlen),info.pos2=as.numeric (info.pos2)) %>%
		mutate (id_sample=paste (id,sample,sep="__")) %>%
		data.frame

	deletions <- delly_data_aux %>%
		filter (info.svtype=="DEL") %>%
		mutate (chrom1=chrom,pos1=pos,chrom2=chrom,pos2=info.end,SVtype=info.svtype,strand1="+",strand2="-") %>%
		dplyr::select (sample,chrom1,pos1,chrom2,pos2,SVtype,strand1,strand2,gt) %>%
		data.frame

	duplications <- delly_data_aux %>%
		filter (info.svtype=="DUP") %>%
		mutate (chrom1=chrom,pos1=pos,chrom2=chrom,pos2=info.end,SVtype=info.svtype,strand1="-",strand2="+") %>%
		dplyr::select (sample,chrom1,pos1,chrom2,pos2,SVtype,strand1,strand2,gt) %>%
		data.frame

	inversion_head_to_head <- delly_data_aux %>%
		filter (info.svtype=="BND",stringr::str_detect (alt,"^\\]|^\\[")) %>%
		mutate (chrom1=chrom,pos1=pos,chrom2=info.chr2,pos2=info.pos2,SVtype="h2hINV",strand1="+",strand2="+") %>%
		dplyr::select (sample,chrom1,pos1,chrom2,pos2,SVtype,strand1,strand2,gt) %>%
		data.frame

	inversion_tail_to_tail <- delly_data_aux %>%
		filter (info.svtype=="BND",stringr::str_detect (alt,"\\]$|\\[$")) %>%
		mutate (chrom1=chrom,pos1=pos,chrom2=info.chr2,pos2=info.pos2,SVtype="t2tINV",strand1="-",strand2="-") %>%
		dplyr::select (sample,chrom1,pos1,chrom2,pos2,SVtype,strand1,strand2,gt) %>%
		data.frame

	result <- rbind (deletions,duplications,inversion_head_to_head,inversion_tail_to_tail)

	return (result)
}

sv_find_event <- function (data, group)
{
	#print ("sv_find_event:")
	#print (data)
	#print (group)

	flank <- 50

	pos1r <- IRanges (start=data %>% pull (pos1),width=1)#IRanges auto merges unless width/end is specified
	pos1rem <- IRanges::reduce (IRanges::resize (pos1r,width=flank*2,fix="center"))

	pos2r <- IRanges (start=data %>% pull (pos2),width=1)#IRanges auto merges unless width/end is specified
	pos2rem <- IRanges::reduce (IRanges::resize (pos2r,width=flank*2,fix="center"))

	# findOverlaps (query, subject
	pos1r_in_pos1rem <- IRanges::findOverlaps (pos1r,pos1rem,type="within",select="all")
	pos2r_in_pos2rem <- IRanges::findOverlaps (pos2r,pos2rem,type="within",select="all")

	result <- data %>%
		mutate (region1=!!S4Vectors::subjectHits (pos1r_in_pos1rem)) %>%
		mutate (region2=!!S4Vectors::subjectHits (pos2r_in_pos2rem)) %>%
		#bind_cols (event_id=group_indices (.,region1,region2)) %>%
		group_by (region1,region2) %>%
		mutate (event_id=cur_group_id ()) %>%
		data.frame

	return (result)
}

as_iranges <- function (data)
{
        mdata <- data %>% dplyr::select (-c (start,end,width)) %>% data.frame
        result <- IRanges (start=data$start,end=data$end,width=data$width)
        mcols (result) <- mdata
        return (result)
}

as_tibble <- function (data)
{
	result <- tibble (
		start = start (data),
		end = end (data),
		width = width (data)
	)
	if ( "revmap" %in% names (mcols (data)) )
	{
		attr(result,"revmap") <- mcols (data)$revmap
	}
	return (result)
}

cleanup_cnv_overlaps <- function (data,group)
{
	#print ("cleanup_cnvs")
	#print (head(data))
	#print (group)

	result <- data %>%
		mutate (width=end-start+1) %>%
		as_iranges %>%
		IRanges::reduce () %>%
		as_tibble %>%
		data.frame

	return (result)
}

cleanup_cnv_disjoin <- function (data, group)
{
	#print (data)
	#print (nrow(data))
	data_with_id <- data %>%
		mutate (id=row_number ())

	data_disjoin <- data_with_id %>%
		as_iranges %>%
		IRanges::disjoin (with.revmap=T) %>%
		as_tibble

	#print ("data_disjoin")
	#print (data_disjoin)
	# Annotate meta data using revmap
	disjoin_id_map <- extractList (data_with_id %>% pull (id), data_disjoin %>% attr("revmap") )

	#print ("disjoin_id_map:")
	#print (disjoin_id_map)

	data_disjoin_id_map <- do.call (rbind,lapply (seq_along (disjoin_id_map), function (x,a) { data.frame (id=x,value=a[[x]]) },a=disjoin_id_map))
	#print (data_disjoin_id_map)

	result <- data_disjoin %>%
		mutate(id=row_number ()) %>%
		inner_join (data_disjoin_id_map,by=c("id"="id")) %>%
		inner_join (data_with_id,by=c("value"="id"),suffix=c("",".original")) %>%
		filter (width>1) %>%
		dplyr::select (start,end,cn) %>%
		data.frame

	return (result)
}

cn_from_cnv_kit <- function (data)
{
	result <- data %>%
		mutate (cn=round (2^log2,digits=0)) %>%
		dplyr::select (sample,chrom,start,end,cn) %>%
		group_by (sample,chrom,cn) %>%
		group_modify (cleanup_cnv_overlaps) %>%
		ungroup () %>%
		group_by (sample,chrom) %>%
		group_modify (cleanup_cnv_disjoin) %>%
		ungroup () %>%
		arrange (sample,chrom,start) %>%
		data.frame

	return (result)
}

