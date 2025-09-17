library(GenomicFeatures)
library(org.Mm.eg.db)
#library(org.Hs.eg.db)

txdb <- makeTxDbFromGFF("Mus_musculus.GRCm39.110.gtf",format="gtf")

#通过exonsBy获取每个gene上的所有外显子的起始位点和终止位点,用reduce去除掉重叠冗余的部分,最后计算长度
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_len <- lapply(exons_gene,function(x){sum(width(reduce(x)))})

exons_gene_len_m <- as.matrix(exons_gene_len)
exons_gene_len_df <- as.data.frame(exons_gene_len_m)
colnames(exons_gene_len_df) <- c('eff_length')
write.csv(exons_gene_len_m, "mm114_gene_length.csv", row.names = TRUE)
