library("pasilla")
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
counts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]


cts <- as.matrix(read.csv("example_data.csv", header = TRUE, row.names = 1))
condition <- c(rep("treated", 3), rep("untreated", 4))
col <- data.frame(condition)
rownames(col) <- colnames(cts)
cts <- cts[, rownames(col)]
dds <- DESeqDataSetFromMatrix(countData = cts, colData = col, design = ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds$condition <- relevel(dds$condition, "untreated")
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
resSig <- subset(resOrdered, padj < 0.05)
write.csv(as.data.frame(res), file = "DE_results.csv")
write.csv(as.data.frame(resOrdered), file = "DE_ordered_results.csv")
write.csv(as.data.frame(resSig), file = "DE_sig_results.csv")
pdf("MA_plot.pdf")
plotMA(res)
dev.off
pdf("MA_plot_resSig.pdf")
plotMA(resSig)
dev.off