library(tidyr)
library(dplyr)
library(pheatmap)
library(Rtsne)
library(data.table)
library(ggplot2)
library(matrixStats)

Files=list.files(pattern="Positions.txt")
#Files=Files[42:56]

results=list()
for (file in Files)
{
	ID=paste(unlist(strsplit(file,"\\."))[1],unlist(strsplit(file,"\\."))[2],sep="_")
	print(paste0("Getting SNP for ",ID))
	tab=read.delim(file)
	colnames(tab)=c("CHROM","POS","REF","ALT","AF","AD0","AD1", "MMQ","MBQ")
	tab=tbl_df(tab)
	tab=tab %>% 
	filter(AD0+AD1>=10) %>%
	filter(MMQ>=50) %>%
	mutate(REF=as.character(REF)) %>%
	mutate(ALT=as.character(ALT)) %>%
	mutate(CHROM=as.character(CHROM))
	tab$ID=ID
	vec=grep(",",tab$ALT)
	if (length(vec) > 1){
		tab=tab[-grep(",",tab$ALT),]
		vec=grep("Phix",tab$CHROM)
	}
		if (length(vec) > 1){
		tab=tab[-grep("Phix",tab$CHROM),]
	}
	tab$uniquepos=paste(tab$CHROM,tab$POS,tab$REF,tab$ALT,sep="_")
	tab = tab[sample(nrow(tab),6000000,replace=T),]
	tab=unique(tab)
	results[[ID]]=tab
}
results = bind_rows(results)

results = results[,c("uniquepos","AF","ID")]
spreads = spread(results,"ID","AF")

#spreads[is.na(spreads)]=0
df = as.data.frame(spreads)
df = df[sample(nrow(df),1000000),]
rownames(df) = df$uniquepos
df=df[,c(2:ncol(df))]
for (i in 1:ncol(df)) {df[is.na(df[,i]),i]=0}
#spreads=spreads[,c(1,2,4,5,6,8)]
for (i in 1:ncol(df)) {df[,i]=as.numeric(df[,i])}
df=df[rowSums(df)>0,]
df=df[,colSums(df)>0]
#write.table(spreads,"Spreads.txt",row.names=T,col.names=T,quote=F,sep="\t")
saveRDS(df,file="df.RDS")
spreads=readRDS("df.RDS")
vec=colnames(spreads)

pheatmap(cor(spreads, use = "complete.obs"),fontsize_row = 3, fontsize_col = 3)

graphics.off()

pdf("Sample_to_Sample_Correlation.pdf", width=20, height = 20)
pheatmap(cor(spreads, use = "complete.obs"),fontsize_row = 6, fontsize_col = 6)
dev.off()

col=read.table("../SNP_Normal/Annotation.txt",header=T,stringsAsFactors=F)

colours=c("black","green","red","blue","brown","gold","darkblue","darkmagenta")
names(colours)=unique(col[,"Group"])

pdf("Sample_to_Sample_Correlation.pdf")
pheatmap(cor(spreads, use = "complete.obs"),fontsize_row = 3, fontsize_col = 3, annotation_col=col,annotation_row=col,annotation_colors=list(Group=colours))
dev.off()

temp=spreads
temp[temp>0]=1

pdf("Sample_to_Sample_Correlation_normed.pdf")
pheatmap(cor(temp, use = "complete.obs"),fontsize_row = 3, fontsize_col = 3,annotation_col=col,annotation_row=col,annotation_colors=list(Group=colours))
dev.off()

temp=spreads
temp[temp<=0.3]=0
temp[temp>=0.7]=0
temp[temp>0.3 & temp < 0.7]=1
temp=temp[rowSums(temp)>0,]

pdf("Sample_to_Sample_Correlation_hetero.pdf")
pheatmap(cor(temp, use = "complete.obs"),fontsize_row = 3, fontsize_col = 3, annotation_col=col,annotation_row=col,annotation_colors=list(Group=colours))
dev.off()

temp=spreads
temp[temp<=0.7]=0
temp[temp>=0.7]=1
temp=temp[rowSums(temp)>0,]

pdf("Sample_to_Sample_Correlation_homo.pdf")
pheatmap(cor(temp, use = "complete.obs"),fontsize_row = 3, fontsize_col = 3,annotation_col=col,annotation_row=col,annotation_colors=list(Group=colours))
dev.off()

col[,"colours"] <- colours[match(col$Group, names(colours))]

annotation=data.table(col)
annotation[,V1:=rownames(col)]
dt=data.table(spreads)
#annotation <- fread("../../Downloads/Annotation.txt", header=T)
#dt <- readRDS("../../Downloads/Spreads.RDS")

dt[, variance := rowVars(as.matrix(.SD))] # get variance
setorder(dt, -variance)
varDT <- dt[, -c("variance")] # select top x most variant rows (without variance column)

pca <- prcomp(t(varDT))

# get explained variance
percentVar <- round(data.table(summary(pca)$importance)[2,1:2]*100)

# get principal components 1 and 2
pcDT <- as.data.table(pca$x)[,1:2]
pcDT[, sample := rownames(pca$x)]

# get groups
pcDT <- merge(pcDT, annotation, by.x="sample", by.y="V1")

# set colors
group.colors <- c(DS2 = "black", HP = "green", BLMM ="red", C = "blue", DS = "brown",
                  PF = "gold", VBC="darkblue", VPF="darkmagenta")

ggplot(pcDT, aes(PC1, PC2, color=Group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[,1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[,2], "% variance")) +
  theme_classic(base_size = 12) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0, color="grey") +
  scale_color_manual(values=group.colors)

ggsave("PCA_all.pdf", plot = last_plot(), device = "pdf", width = 10, height = 10)

spreads=t(spreads)

set.seed(1)  # set random seed

perplexity=1

# plot 2D t-SNE projection
for (perplexity in 1:50)
{
	pdf(paste(perplexity,".PCA.pdf",sep=""))
	rtsne_out <- Rtsne(spreads, pca = TRUE, verbose = TRUE, dims=2, perplexity=perplexity, check_duplicates=F, num_threads=8, max_iter = 2000)
	plot(rtsne_out$Y, asp = 1, pch = 20,col=col[,"colours"], 
	     cex = 0.75, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5, 
	     xlab = "t-SNE dimension 1", ylab = "t-SNE dimension 2", 
	     main = "2D t-SNE projection")
	#text(rtsne_out$Y,labels=unlist(strsplit(rownames(spreads),"_NormalOnly")),pos=3,cex= 0.7)
	dev.off()
}


# pdf("True_Values_normed.pdf")
# pheatmap(spreads[runif(10000,1,nrow(temp)),],show_rownames=F)
# dev.off()


# pdf("True_Values.pdf")
# pheatmap(spreads[runif(10000,1,nrow(spreads)),],show_rownames=F)
# dev.off()
