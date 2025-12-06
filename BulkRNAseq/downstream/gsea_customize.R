library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(stringr)library(ggplot2)

rm(list = ls())
options(stringsAsFactors = F)
setwd("yourPath")

deseq_data <- read.table("deseq2.tsv", header=TRUE, sep="\t")
colnames(deseq_data)[1] <- c("geneid")
analyse_data <- deseq_data %>% dplyr::select(geneid, log2FoldChange) %>% arrange(desc(log2FoldChange)) %>% as.data.frame

gene_list <- analyse_data$log2FoldChange
names(gene_list) <- analyse_data$geneid
gmt <- read.gmt("yourgmt")
gsea_res <- GSEA(
  geneList = gene_list,
  TERM2GENE = gmt,
  pvalueCutoff = 1,
  #nPermSimple = 10000,
  maxGSSize = 2000
)
summary(gsea_res@result$pvalue)
pv = signif(gsea_res@result$pvalue, 3)
gsea_pic1 <- gseaplot2(gsea_res, title = paste0(gsea_res$Description[1],', p =',pv), geneSetID = 1)
gsea_pic1
gsea_pic2 <- gseaplot(gsea_res, title = paste0(gsea_res$Description[1],', p=',pv), geneSetID = 1, by = "runningScore", pvalue_table = T) +
  theme(panel.grid=element_blank())
gsea_pic2
ggsave(filename = paste0('gsea_',gsea_res$Description[1],'_2.pdf'), plot=gsea_pic1, width = 8, height = 6, units = "in", dpi = 300)
ggsave(filename = paste0('gsea_',gsea_res$Description[1],'_1.pdf'), plot=gsea_pic2, width = 8, height = 4.5, units = "in", dpi = 300)
