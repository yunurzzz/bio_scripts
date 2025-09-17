library(DESeq2)
library(ggplot2)
library(ggthemes)
library(ggrepel)

rm(list = ls())
options(stringsAsFactors = F)
setwd("yourPath")

dat <- read.table('count.tsv', sep = '\t', row.names = 1, header = 1, check.names = FALSE)

condition <- factor(c(rep("treat",3),rep("ctrl",3)))
coldata <- data.frame(row.names=colnames(dat), condition)
dds1 <- DESeqDataSetFromMatrix(countData = dat,
                               colData = coldata,
                               design = ~condition)


###PCA
vsd <- vst(dds1, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c('condition'), returnData = T)
plotPCA(vsd, intgroup=c('condition'), returnData = F)
pic1 <- plotPCA(vsd, intgroup=c('condition'), returnData = F) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
pic1

pic2 <- ggplot(pcaData, aes(x=PC1, y=PC2), color = condition, shape = group)+
  geom_point(aes(color = condition), size = 3.5) +
  scale_color_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf')) +
  xlab("PC1: %") + ylab("PC2: %") +
  theme_bw() + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  geom_text_repel(aes(label = name), size = 2.8, show.legend = FALSE,
                  box.padding = unit(0.5, 'lines'))
pic2

ggsave(filename = "./pca.pdf", plot=pic2, width = 6, height = 4.2, units = "in", dpi = 300)

