library(DESeq2)
library(ggplot2)
library(ggthemes)
library(org.Mm.eg.db)
#library(org.Hs.eg.db)
library(clusterProfiler)

rm(list = ls())
options(stringsAsFactors = F)
setwd("yourPath")

dat <- read.table('count.tsv', sep = '\t', row.names = 1, header = 1, check.names = FALSE)
nrow(dat)
dat <- dat[rowSums(dat)>0,]
nrow(dat)
dat_adj <- dat[,1:ncol(dat)]+1

condition <- factor(c(rep("treat",3),rep("ctrl",3)))

coldata <- data.frame(row.names=colnames(dat_adj), condition)
dds1 <- DESeqDataSetFromMatrix(countData = dat_adj,
                               colData = coldata,
                               design = ~condition)
vsd <- vst(dds1) 
plotPCA(vsd, intgroup=c('condition')) + theme_bw()

dds2 <- DESeq(dds1)
#注意，需将treat在前，control在后，意为treat相较于control中哪些基因上调/下调
res <- results(dds2, contrast = c('condition', 'treat', 'ctrl'))
summary(res)
res_out <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

##筛选差异表达基因
#首先对表格排个序，按padj值升序排序，相同padj值下继续按log2FC降序排序
res_out <- res_out[order(res_out$padj, res_out$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
write.table(res_out, './deseq2.tsv', col.names = NA, sep = '\t', quote = FALSE)

#log2FC≥1&padj<0.01标识up，代表显著上调的基因
#log2FC≤-1&padj<0.01标识down，代表显著下调的基因
#其余标识 none，代表非差异的基因
res_out[which(res_out$log2FoldChange > 0 & res_out$padj < 0.05),'sig'] <- 'up'
res_out[which(res_out$log2FoldChange < 0 & res_out$padj < 0.05),'sig'] <- 'down'
res_out[which(res_out$padj >= 0.05),'sig'] <- 'none'
res_out_select <- subset(res_out, sig %in% c('up', 'down'))
#根据up和down分开输出
res_out_up <- subset(res_out, sig == 'up')
res_out_down <- subset(res_out, sig == 'down')
# 防止输出错位
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}
res_out_up <- adjustdata (res_out_up)
res_out_down <- adjustdata (res_out_down)
write.table(res_out_up, file = './deseq2_up.tsv', sep = '\t', row.names = F, col.names = T, quote = FALSE)
write.table(res_out_down, file = './deseq2_down.tsv', sep = '\t', row.names = F, col.names = T, quote = FALSE)

# geneid2symbol
res_out2 = res_out
res_out$ENSEMBL = rownames(res_out)
geneid = res_out$ENSEMBL
transid <- bitr(geneid, OrgDb = org.Mm.eg.db, toType = c("SYMBOL","GENENAME"), fromType = "ENSEMBL")
res_convert <- merge(res_out,transid,by='ENSEMBL')
col_n <- ncol(res_convert)
res_convert <- res_convert[,c((col_n-1),2:(col_n-2),col_n)]
res_convert <- res_convert[order(res_convert$padj, res_convert$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
write.table(res_convert, file = './deseq2_symbol.tsv', sep = '\t', row.names = F, col.names = T, quote = FALSE)
res_convert[which(res_convert$log2FoldChange > 0 & res_convert$padj < 0.05),'sig'] <- 'up'
res_convert[which(res_convert$log2FoldChange < 0 & res_convert$padj < 0.05),'sig'] <- 'down'
res_convert[which(res_convert$padj >= 0.05),'sig'] <- 'none'
col_n <- ncol(res_convert)
res_convert_up <- subset(res_convert, sig == 'up')
res_convert_down <- subset(res_convert, sig == 'down')
write.table(res_convert_up, file = './deseq2_symbol_up.tsv', sep = '\t', row.names = F, col.names = T, quote = FALSE)
write.table(res_convert_down, file = './deseq2_symbol_down.tsv', sep = '\t', row.names = F, col.names = T, quote = FALSE)

