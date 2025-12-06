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

###gseGO
gsego_res <- gseGO(
  geneList = gene_list,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "ALL",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  eps = 1e-10,
  verbose = FALSE
  )

gsego_res_all = as.data.frame(gsego_res@result)
nrow(gsego_res_all)
gene_ratio = as.data.frame(str_count(gsego_res_all$core_enrichment,pattern = "/"))
colnames(gene_ratio) <- c("genecount")
gsego_res_all$geneCount = gene_ratio$genecount
gene_ratio$gene_ratio = gene_ratio[,1]/gsego_res_all[,4]
gsego_res_all$geneRatio = gene_ratio$gene_ratio
gsego_result_sorted = gsego_res_all[order(gsego_res_all$qvalue),]
gsego_result_sorted$qvalue = signif(gsego_result_sorted$qvalue,3)

gsego_result_sorted[which(gsego_result_sorted$NES > 0),'sig'] <- 'up'
gsego_result_sorted[which(gsego_result_sorted$NES < 0),'sig'] <- 'down'
gsego_up <- subset(gsego_result_sorted, sig == 'up')
gsego_down <- subset(gsego_result_sorted, sig == 'down')
write.table(gsego_up, file = 'gesgo_up.tsv', sep = '\t', row.names = F, quote = FALSE)
write.table(gsego_down, file = 'gesgo_down.tsv', sep = '\t', row.names = F, quote = FALSE)
