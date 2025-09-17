library(ggplot2)
library(ggrepel)
library(DESeq2)
library(limma)

rm(list = ls())
options(stringsAsFactors = F)
setwd("yourPath")

dat <- read.delim('count.txt', row.names = 1, sep = '\t', check.names = FALSE)
coldata <- read.delim('group.txt', row.names = 1, sep = '\t', check.names = FALSE)
dds1 <- DESeqDataSetFromMatrix(countData = dat,
                               colData = coldata,
                               design = ~batch+condition)

vsd <- vst(dds1, blind=FALSE)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
pcaData <- plotPCA(vsd, intgroup=c("condition", "batch"), returnData=TRUE)
plotPCA(vsd, intgroup=c("condition", "batch"), returnData=F)
pcaData <- plotPCA(vsd, intgroup=c("condition", "batch"), returnData=T)
write.table(pcaData, 'pca_data.tsv', col.names = NA, sep = '\t', quote = FALSE)

pcaData = read.table('pca_data.tsv',header = T, row.names = NULL, check.names = FALSE, sep = '\t')
pic2 <- ggplot(pcaData, aes(PC1, PC2, color=condition)) + 
  geom_point(aes(color = condition), size = 3) +
  scale_color_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf')) +
  #scale_color_manual(values = c('#1b9e77', '#d95f02')) +
  xlab("PC1: 57%") + ylab("PC2: 19%") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())#+
  geom_text_repel(aes(label = name), size = 3, show.legend = FALSE,
                  box.padding = unit(0.5, 'lines'), max.overlaps = Inf)
pic2
ggsave(filename = "./pca.pdf", plot=pic2, width = 6, height = 3.8, units = "in", dpi = 300)
