library(stringr)
library(dplyr)
library(org.Mm.eg.db)
#library(org.Hs.eg.db)
library(clusterProfiler)
#library(xlsx)

rm(list = ls())
options(stringsAsFactors = F)
setwd("yourPath")

#要保证表达矩阵的行名和存放基因长度向量的名字一致, 这一步非常重要
expMatrix = read.table('count.tsv', row.names = 1, header = T, sep = '\t', check.names=F)
expMatrix[1:3,]
eff_length <- read.csv("mm110_gene_length.csv", header = T)
rownames(eff_length) <- eff_length$gene_id
colnames(eff_length) <- c("gene_id","eff_length")
rownames(eff_length) <- do.call(rbind,strsplit(as.character(eff_length$gene_id),'\\.'))[,1]

feature_ids <- rownames(expMatrix)

#检查gtf eff_length文件和表达量输入文件里基因名的一致性
if (! all(feature_ids %in% rownames(eff_length))){
  tbl <- table(feature_ids %in% rownames(eff_length))
  msg1 <- sprintf("%i gene is shared, %i gene is specified", tbl[[2]],tbl[[1]])
  warning(msg1)
}

if (! identical(feature_ids, rownames(eff_length))){
  msg2 <- sprintf("Given GTF file only contain %i gene, but experssion matrix has %i gene", nrow(eff_length), nrow(expMatrix))
  warning(msg2)
}

# trim the expression matrix and effetive gene length
expMatrix <- expMatrix[feature_ids %in% rownames(eff_length),]
mm <- match(rownames(expMatrix), rownames(eff_length))
eff_length <- eff_length[mm, ]

if (identical(rownames(eff_length), rownames(expMatrix))){
  print("GTF and expression matix now have the same gene and gene in same order")
}

countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp(log(counts) + log(1e9) - log(effLen) - log(N))
}
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

#FPKM
fpkms <- apply(expMatrix, 2, countToFpkm, effLen = eff_length$eff_length)
fpkms.m <- data.frame(round(fpkms,3))
colnames(fpkms.m) <- colnames(expMatrix)
dim(fpkms.m)
fpkms.m[1:3,]
write.table(fpkms.m, file="fpkm.tsv", sep="\t", quote=F, row.names=T)
fpkms.m$ENSEMBL <- rownames(fpkms.m)
geneid = fpkms.m$ENSEMBL
transid <- bitr(geneid, OrgDb = org.Mm.eg.db, toType = c("SYMBOL"), fromType = "ENSEMBL")
hitmap <- merge(fpkms.m,transid,by='ENSEMBL')
nrow(hitmap)
nrow(fpkms.m)
nrow(transid)
col_n <- ncol(hitmap)
hitmap <- hitmap[,c(col_n,2:(col_n-1))]
write.table(hitmap, file="fpkm_symbol.tsv", sep="\t", quote=F, row.names=F)

#####TPM
tpms <- apply(expMatrix, 2, countToTpm, effLen = eff_length$eff_length)
tpms.m <- data.frame(round(tpms,3))
colnames(tpms.m) <- colnames(expMatrix)
dim(tpms.m)
tpms.m[1:3,]
write.table(tpms.m, file="tpm.tsv", sep="\t", quote=F, row.names=T)
tpms.m$ENSEMBL <- rownames(tpms.m)
geneid = tpms.m$ENSEMBL
transid <- bitr(geneid, OrgDb = org.Mm.eg.db, toType = c("SYMBOL"), fromType = "ENSEMBL")
hitmap <- merge(tpms.m,transid,by='ENSEMBL')
nrow(hitmap)
nrow(tpms.m)
nrow(transid)
col_n <- ncol(hitmap)
hitmap <- hitmap[,c(col_n,2:(col_n-1))]
write.table(hitmap, file="tpm_symbol.tsv", sep="\t", quote=F, row.names=F)
