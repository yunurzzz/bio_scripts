library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(stringr)

rm(list = ls())
options(stringsAsFactors = F)
setwd("yourPath")

deseq_data <- read.table("deseq2.tsv",header=TRUE,sep="\t")
colnames(deseq_data)[1] <- c("geneid")
analyse_data <- deseq_data %>% dplyr::select(geneid, log2FoldChange) %>% arrange(desc(log2FoldChange)) %>% as.data.frame
gene_list <- analyse_data$log2FoldChange
names(gene_list) <- analyse_data$geneid

###gseKEGG
transid <- bitr(analyse_data$geneid, OrgDb = org.Mm.eg.db, toType = "ENTREZID", fromType = "ENSEMBL")
analyse_data <- analyse_data %>% inner_join(transid,by=c("geneid"="ENSEMBL"))
analyse_data <- analyse_data %>% arrange(desc(log2FoldChange),multiple = "all")
gene_list <- analyse_data$log2FoldChange
names(gene_list) <- analyse_data$ENTREZID
gsekegg_res <- gseKEGG(
  geneList = gene_list,
  organism = 'mmu',
  #keyType = "ENTREZID",
  pvalueCutoff = 1,
  eps = 1e-10,
  verbose = F
  )
gsekegg_res_all = as.data.frame(gsekegg_res@result)
nrow(gsekegg_res_all)
gene_ratio = as.data.frame(str_count(gsekegg_res_all$core_enrichment,pattern = "/"))
colnames(gene_ratio) <- c("genecount")
gsekegg_res_all$geneCount = gene_ratio$genecount
gene_ratio$gene_ratio = gene_ratio[,1]/gsekegg_res_all[,3]
gsekegg_res_all$geneRatio = gene_ratio$gene_ratio
gsekegg_result_sorted = gsekegg_res_all[order(gsekegg_res_all$qvalue),]
gsekegg_result_sorted$qvalue = signif(gsekegg_result_sorted$qvalue,3)
gsekegg_result_sorted$Description = sub("-[^-]+$","",gsekegg_result_sorted$Description)

gsekegg_result_sorted[which(gsekegg_result_sorted$NES > 0),'sig'] <- 'up'
gsekegg_result_sorted[which(gsekegg_result_sorted$NES < 0),'sig'] <- 'down'
gsekegg_up <- subset(gsekegg_result_sorted, sig == 'up')
gsekegg_down <- subset(gsekegg_result_sorted, sig == 'down')
write.table(gsekegg_up, file = 'gsekegg_up.tsv', sep = '\t', row.names = F, quote = FALSE)
write.table(gsekegg_down, file = 'gsekegg_down.tsv', sep = '\t', row.names = F, quote = FALSE)
