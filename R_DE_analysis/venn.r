install.packages("gplots")
library(gplots)
DESeq2_top10 <- resSig@rownames[1:10]
NOISeq_top10 <- row.names(results.deg)[1:10]
t <-as.data.frame(topTags(lrt))
edgeR_top10 <- row.names(t)
pdf("DE_Venn.pdf")
venn(list(DESeq2 = DESeq2_top10, edgeR = edgeR_top10, NOIseq = NOISeq_top10))
dev.off()
all < - Reduce(intersect, list(DESeq2_top10, NOISeq_top10, edgeR_top10))
d_vs_n < - Reduce(intersect, list(DESeq2_top10, NOISeq_top10))
d_vs_e < - Reduce(intersect, list(DESeq2_top10, edgeR_top10))
n_vs_e < - Reduce(intersect, list(NOISeq_top10, edgeR_top10))
