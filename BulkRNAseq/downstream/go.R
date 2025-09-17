library(stringr)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
#library(org.Hs.eg.db)
library(forcats)
library(enrichplot)

rm(list = ls())
options(stringsAsFactors = F)
setwd("yourPath")

rt = read.table('deseq2_up.tsv', header = T, sep = '\t')
colnames(rt)[1] <- c("geneid")
geneid = rt$geneid
head(geneid)

transid <- bitr(geneid, OrgDb = org.Mm.eg.db, toType = c("GENENAME","SYMBOL","ENTREZID"), fromType = "ENSEMBL")

ego <- enrichGO(gene = transid$ENTREZID, #ENSEMBL
                OrgDb = org.Mm.eg.db,
                #OrgDb = org.Hs.eg.db,
                keyType = 'ENTREZID', #ENSEMBL
                ont= "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                readable = T)

# Get the similarity matrix
ego1 <- pairwise_termsim(ego)
emapplot(ego1)
ego2 <- simplify(ego1, cutoff=0.7, by="p.adjust", select_fun=min)

#获取所有分析结果
#go_result_all = ego1@result
go_result_all = ego2@result
go_result = as.data.frame(go_result_all)

#将GENERATIO结果转化为小数
gene_ratio = as.data.frame(apply(str_split(go_result$GeneRatio,"/",simplify=T),2,as.numeric))
gene_ratio$gene_ratio = gene_ratio[,1]/gene_ratio[,2]
go_result$GeneRatio = gene_ratio$gene_ratio

#按照qvalue排序
go_result_sorted = go_result[order(go_result$qvalue),]

#结果写入表格
write.table(go_result_sorted,file="go_up.tsv",sep = "\t",quote = F,row.names = F)
