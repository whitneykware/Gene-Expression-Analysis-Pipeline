#DESeq2
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library('DESeq2')
setwd("~/RNASeq/ReadCounts")
directory = '/Users/whitneyware/RNASeq/ReadCounts'

#from HTSeq files
sampleFiles <- grep('counts.txt', list.files(directory), value=TRUE)
sampleNames <- sub('_counts.txt','',sampleFiles)
sampleNames <- sub('.*\\_', '', sampleNames)
sampleCondition <- c(rep("SPF",9),rep("GF",7))
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, ~ condition)

#from count matrix (unnormalized)
cts<-as.matrix(read.csv("MyData.csv", header = TRUE, row.names = 1))
condition<-c(rep("treated",3),rep("untreated",4))
coldata<-data.frame(condition)
rownames(coldata)<-colnames(cts)
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)

#filter low counts
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds$condition <- factor(dds$condition, levels = c("untreated","treated"))

# Differential expression analysis
dds <- DESeq(dds)

#results
res <- results(dds)
resOrdered<- res[order(res$padj),]
resSig <- subset(resOrdered, padj<0.05)

#write results to file
write.csv(as.data.frame(res), file = "DESeq2_DE_Results.csv")
write.csv(as.data.frame(resOrdered), file = "DESeq2_DE_Ordered_Results.csv")
write.csv(as.data.frame(resSig), file = "DESeq2_DE_Sig_Results.csv")

#Plot showing differentially expressed genes
pdf("MA_plot.pdf")
plotMA(resSig)
dev.off()

pdf("MA_plot_resSig.pdf")
plotMA(resSig)
dev.off()

