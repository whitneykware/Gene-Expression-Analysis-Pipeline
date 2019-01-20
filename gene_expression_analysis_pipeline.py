
class edge_R:
    def __init__(self, filename):
        self.filename = filename
        self.exp_group_name = input("Enter experimental group name: ")
        self.exp_sample_number = int(input("Enter the number of experimental samples: "))
        self.control_group = input("Enter control group name: ")
        self.con_sample_number = int(input("Enter the number of control samples: "))

    def edge_R_script(self):
        with open(self.filename, 'r') as inf, open('edge_r_out.r', 'w') as outf:
            outf.write("library(edgeR)\n")
            outf.write("x <- read.csv('%s', row.names = 1)\n" % self.filename)
            outf.write("group <- factor(c(rep('%s',%d), rep('%s',%d)))\n" %
                       (self.exp_group_name, self.exp_sample_number, self.control_group, self.con_sample_number))
            outf.write("y <- DGEList(counts = x,group=group)\n")
            outf.write("y <- calcNormFactors(y)\n")
            outf.write("design <- model.matrix(~group)\n")
            outf.write("y<-estimateDisp(y,design)\n")
            outf.write("fit <-glmFit(y,design)\n")
            outf.write("lrt <-glmLRT(fit,coef=2)\n")
            outf.write("topTags(lrt)\n")
            outf.write("write.csv(topTags(lrt),file= 'top_edgeR_lrt.csv')\n")
            outf.write("is.de <- decideTestsDGE(lrt)\n")
            outf.write("summary(is.de)\n")
            outf.write("plotMD(lrt,status = is.de,values = c(1,-1),col=c('red','blue'),legend = 'topright')\n")
        outf.close()


class DESeq2:

    def __init__(self):
        self.datatype = input('Input format: HTSeq files or count matrix? ')
        if self.datatype == "count matrix":
            self.matrix = input('Enter count matrix filename: ')
        else:
            self.files = "counts"
            self.countFilesPath = input('Enter path to count files: ')
        self.controlName = input('Name of Control Group: ')
        self.controlSamples = int(input('Number of Control Samples: '))
        self.treatmentName = input('Name of Experimental Group: ')
        self.treatmentSamples = int(input('Number of Experimental Samples: '))

    def create_deseq2_script(self):
        with open('deseq2_DE_script.r', 'w') as outf:
            outf.write('source("https://bioconductor.org/biocLite.R")\n'+'biocLite("DESeq2")\n'
                    +'library("DESeq2")\n')
            if self.datatype == "HTSeq files":
                outf.write('directory = "{}"\n'.format(self.countFilesPath))
                outf.write('sampleFiles <- grep("{}.txt", list.files(directory), value = TRUE)\n'.format(self.files))
                outf.write('sampleNames <- sub("_counts.txt", "", sampleFiles)\n')
                outf.write('sampleNames <- sub(".*\\_", "", sampleNames)\n')
                outf.write('sampleCondition <- c(rep("{0}", {1}), rep("{2}", {3}))\n'.format(self.treatmentName,
                            self.treatmentSamples, self.controlName, self.controlSamples))
                outf.write('sampleTable <- data.frame(sampleName = sampleNames, '
                       'fileName = sampleFiles, condition = sampleCondition)\n')
                outf.write('dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, '
                       'directory = directory, design = ~ condition)\n')
            elif self.datatype == "count matrix":
                outf.write('cts <- as.matrix(read.csv("{}", header = TRUE, '
                           'row.names = 1))\n'.format(self.matrix))
                outf.write('condition <- c(rep("{0}", {1}), rep("{2}", {3}))\n'.format(self.treatmentName,
                            self.treatmentSamples, self.controlName, self.controlSamples))
                outf.write('col <- data.frame(condition)\n')
                outf.write('rownames(col) <- colnames(cts)\n')
                outf.write('cts <- cts[, rownames(col)]\n')
                outf.write('dds <- DESeqDataSetFromMatrix(countData = cts, '
                           'colData = col, design = ~ condition)\n')
            outf.write('dds <- dds[ rowSums(counts(dds)) > 1, ]\n')
            outf.write('dds$condition <- relevel(dds$condition, "{}")\n'.format(self.controlName))
            outf.write('dds <- DESeq(dds)\n')
            outf.write('res <- results(dds)\n')
            outf.write('resOrdered <- res[order(res$padj),]\n')
            outf.write('resSig <- subset(resOrdered, padj < 0.05)\n')
            outf.write('write.csv(as.data.frame(res), file = "DESeq2_DE_Results.csv")\n')
            outf.write('write.csv(as.data.frame(resOrdered), file = "DESeq2_DE_Ordered_Results.csv")\n')
            outf.write('write.csv(as.data.frame(resSig), file = "DESeq2_DE_Sig_Results.csv")\n')
            outf.write('pdf("MA_plot.pdf")\n')
            outf.write('plotMA(res)\n')
            outf.write('dev.off()\n')
            outf.write('pdf("MA_plot_resSig.pdf")\n')
            outf.write('plotMA(resSig)\n')
            outf.write('dev.off()\n')
        outf.close()


class NOISeq:

    def __init__(self):
        self.matrix = input('Enter count matrix filename: ')
        self.controlName = input('Name of Control Group: ')
        self.controlSamples = int(input('Number of Control Samples: '))
        self.treatmentName = input('Name of Experimental Group: ')
        self.treatmentSamples = int(input('Number of Experimental Samples: '))
        self.replicates = input('Replicates?: technical, biological, or no: ')
        self.norm = input('If the data is pre normalized enter n, '
                          'otherwise enter preferred normalization method: (rpkm, uqua, or tmm) ')

    def create_noiseq_script(self):
        with open('noiseq_DE_script.r', 'w') as outf:
            outf.write('source("https://bioconductor.org/biocLite.R")\n'+'biocLite("NOISeq")\n'+'library("NOISeq")\n')
            outf.write('countdata <- as.matrix(read.csv("{}", header = TRUE, row.names = 1))\n'.format(self.matrix))
            outf.write('condition <- c(rep("{0}", {1}), rep("{2}", {3}))\n'.format(self.treatmentName,
                        self.treatmentSamples, self.controlName, self.controlSamples))
            outf.write('coldata <- data.frame(condition)\n')
            outf.write('rownames(coldata) <- colnames(countdata)\n')
            outf.write('countdata<- filtered.data(countdata, factor = coldata$condition, '
                       'norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1, p.adj = "fdr")\n')
            outf.write('mydata <- readData(data = countdata, factors = coldata)\n')
            outf.write('results <- noiseq(mydata, factor = "condition", '
                       'norm = "{0}", replicates = "{1}")\n'.format(self.norm, self.replicates))
            outf.write('results.deg <- degenes(results, q = 0.8, M = NULL)\n')
            outf.write('out <- results@results[[1]]\n')
            outf.write('write.csv(out, file = "NOISeq_DE_Results.csv", quote = FALSE)\n')
            outf.write('write.csv(results.deg, file = "NOISeq_DE_Sig_Results.csv", quote = FALSE)\n')
            outf.write('pdf("EXP_plot.pdf")\n')
            outf.write('DE.plot(results, q=0.9, graphic="expr", log.scale = TRUE)\n')
            outf.write('dev.off()\n')
            outf.write('pdf("MD_plot.pdf")\n')
            outf.write('DE.plot(results, q=0.8, graphic="MD")\n')
            outf.write('dev.off()\n')
        outf.close()


class CompareDEG:

    def __init__(self):
        self.num = int(input('Number of top significant genes to be compared: '))

    def create_venn(self):
        with open('venn.r', 'w') as file:
            file.write('install.packages("gplots")\n')
            file.write('library(gplots)\n')
            file.write('DESeq2_top{0} <- resSig@rownames[1:{1}]\n'.format(self.num, self.num))
            file.write('NOISeq_top{0} <- row.names(results.deg)[1:{1}]\n'.format(self.num, self.num))
            file.write('t <-as.data.frame(topTags(lrt))\n')
            file.write('edgeR_top{} <- row.names(t)\n'.format(self.num))
            file.write('pdf("DE_Venn.pdf")\n')
            file.write('venn(list(DESeq2 = DESeq2_top{0}, edgeR = edgeR_top{1}, NOIseq = NOISeq_top{2}))\n'
                       .format(self.num, self.num, self.num))
            file.write('dev.off()\n')
            file.write('all <- Reduce(intersect, list(DESeq2_top10, NOISeq_top10, edgeR_top10))\n')
            file.write('d_vs_n <- Reduce(intersect, list(DESeq2_top10, NOISeq_top10))\n')
            file.write('d_vs_e <- Reduce(intersect, list(DESeq2_top10, edgeR_top10))\n')
            file.write('n_vs_e <- Reduce(intersect, list(NOISeq_top10, edgeR_top10))\n')
        file.close()


package_choice = input('Would you like to run DESeq2, NOISeq, edgeR, or all: ')

if package_choice == "DESeq2":
    x = DESeq2()
    x.create_deseq2_script()
elif package_choice == "NOISeq":
    x = NOISeq()
    x.create_noiseq_script()
elif package_choice == "edgeR":
    filename = input('Enter filename: ')
    x = edge_R(filename)
    x.edge_R_script()
elif package_choice == "all":
    filename = input('Enter filename: ')
    x = DESeq2()
    x.create_deseq2_script()
    a = NOISeq()
    a.create_noiseq_script()
    b = edge_R(filename)
    b.edge_R_script()
    c = CompareDEG()
    c.create_venn()



