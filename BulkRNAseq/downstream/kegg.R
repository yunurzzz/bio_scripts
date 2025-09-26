library(stringr)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
#library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(enrichplot)

rm(list = ls())
options(stringsAsFactors = F)
setwd("yourPath")

rt = read.table('deseq2_up.tsv',header = T,sep = '\t')
colnames(rt)[1] <- c("geneid")
geneid = rt$geneid
head(geneid)

transid <- bitr(geneid, OrgDb = org.Mm.eg.db, toType = c("ENTREZID","SYMBOL","GENENAME"), fromType = "ENSEMBL")

ekegg <- enrichKEGG(gene = transid$ENTREZID,
                    organism = 'mmu', #hsa
                    #keyType = 'ENTREZID',
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1, #0.01
                    qvalueCutoff = 1) #0.05

ekegg <- setReadable(ekegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")

kegg_result_all = ekegg@result
kegg_result_all$Description=factor(kegg_result_all$Description,levels = kegg_result_all$Description)

kegg_result = as.data.frame(kegg_result_all)

gene_ratio = as.data.frame(apply(str_split(kegg_result$GeneRatio,"/",simplify=T),2,as.numeric))
gene_ratio$gene_ratio = gene_ratio[,1]/gene_ratio[,2]
kegg_result$GeneRatio = gene_ratio$gene_ratio

kegg_result_sorted = kegg_result[order(kegg_result$qvalue),]
kegg_result_sorted$qvaluedigit = signif(kegg_result_sorted$qvalue,3)
kegg_result_sorted$Description = sub("-[^-]*$","",kegg_result_sorted$Description)

#结果写入表格
write.table(kegg_result_sorted,file="kegg_up.tsv",sep = "\t",quote = F,row.names = F)
