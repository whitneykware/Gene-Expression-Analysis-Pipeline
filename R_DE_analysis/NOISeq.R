#NOISeq
source("https://bioconductor.org/biocLite.R")
biocLite("NOISeq")
library(NOISeq)
setwd("~/RNASeq/ReadCounts")

#read in data
countdata<-as.matrix(read.csv("MyData.csv", header = TRUE, row.names = 1))
condition<-c(rep("SPF",9),rep("GF",7))
coldata<-data.frame(condition)
rownames(coldata)<-colnames(countdata)

#make noiseq object
mydata <- readData(data = counts, factors = coldata)

#DE
results <- noiseq(mydata, factor = "condition", norm = "rpkm", replicates = "no")

n<- noiseq(mydata, norm = "rpkm", factor = "condition", replicates = "no")
t <- noiseq(mydata, norm = "tmm", factor = "condition", replicates = "no")
b <- noiseq(mydata, factor = "condition",norm = "n", replicates = "biological")
mynoiseqbio <- noiseqbio(mydata, factor = "condition", adj = 1.5, plot = FALSE, a0per = 0.9, random.seed = 12345, filter = 1)

#results
results.deg <- degenes(results, q = 0.8, M = NULL)
out <- results@results[[1]]

mynoiseqbio.deg <- degenes(mynoiseqbio, q = 0.8, M = NULL)

#write results to file 
write.csv(out, file = "NOISeq_DE_Results.csv", quote = FALSE)
write.csv(results.deg, file = "NOISeq_DE_Sig_Results.csv", quote = FALSE)

#plot
pdf("plot.pdf")
DE.plot(results, q = 0.8, graphic = "expr", log.scale = TRUE)
dev.off()
DE.plot(results, q = 0.8, graphic = "MD")





