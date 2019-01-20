class edge_R:
    def __init__(self, filename):
        self.filename = filename
        self.exp_group_name = input("Enter experimental group name: ")
        self.exp_sample_number = int(input("Enter the number of experimental samples: "))
        self.control_group = input("Enter control group name: ")
        self.con_sample_number = int(input("Enter the number of control samples: "))
        self.norm_factor = input(
            "Enter the preferred method for normalization, from the given options TMM, RLE, upperquartile, none: ")
        self.top_results = int(input("Enter the number of top DE results required: "))
        self.adjust_method = input(
            "Enter the desired PValue adjustment method from BH, fdr, BY, holm, none: ")
        self.sort_top_results = input("Enter if the top results be sorted by PValue, logFC, none: ")
        self.desired_pvalue = float(input("Enter the cutoff value for adjusted pvalues: "))

    def edge_R_script(self):
        with open(self.filename, 'r') as inf, open('edge_r_out.r', 'w') as outf:
            outf.write("source('http://bioconductor.org/biocLite.R')\n")
            outf.write("library(edgeR)\n")
            outf.write("x <- read.csv('%s', row.names = 1)\n" % self.filename)
            outf.write("group <- factor(c(rep('%s',%d), rep('%s',%d)))\n" %
                       (self.exp_group_name, self.exp_sample_number, self.control_group, self.con_sample_number))
            outf.write("y <- DGEList(counts = x,group=group)\n")
            outf.write("y <- calcNormFactors(y, method = '%s')\n" % self.norm_factor)
            outf.write("design <- model.matrix(~group)\n")
            outf.write("y<-estimateDisp(y,design)\n")
            outf.write("fit <-glmFit(y,design)\n")
            outf.write("lrt <-glmLRT(fit,coef=2)\n")
            outf.write("topTags(lrt, n=%d, adjust.method = '%s', sort.by = '%s', p.value = %f)\n" % (
                self.top_results, self.adjust_method, self.sort_top_results, self.desired_pvalue))
            outf.write("write.csv(topTags(lrt),file= 'top_edgeR_lrt.csv')\n")
            outf.write("is.de <- decideTestsDGE(lrt)\n")
            outf.write("summary(is.de)\n")
            outf.write(
                "plotMD(lrt,status = is.de,values = c(1,-1),col=c('red','blue'),legend = 'topright')\n")


filename = input("Enter the filename: ")
e = edge_R(filename)
e.edge_R_script()
