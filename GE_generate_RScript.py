class Pipeline:
    def __init__(self,filename):
        # Taking in filename when creating instance and creating outfile name.
        self.filename = filename
        self.outname = self.filename+"_out.txt"

    def parser(self):
        with open(self.filename, 'r') as infile, open(self.outname,'w') as outfile:
            for line in infile.readlines():
                if line.startswith("!"):
                    continue

                else:
                    outfile.write(line)

    def datainfo(self):
        self.col = []
        self.samplesize = input("Total sample size: ")
        self.group = input("Number of experimental groups: ")
        for i in range(1,int(self.group)+1):
            cols = input("Columns of the %d experimental group:" % i)
            self.col.append(cols)
        for i in range(int(self.group)):
            self.col[i] = self.col[i].split(",")

        self.control = input("Columns of the control group: ")
        self.control = self.control.split(",")
        #print(self.col,self.group,self.control)

    def fillline(self):
        self.fill = [0]*int(self.samplesize)
        for i in range(len(self.col)):
            for j in range(len(self.col[i])):
                self.fill[int(self.col[i][j])-1] = i+1

        for i in range(len(self.control)):
            self.fill[int(self.control[i])-1] = int(self.group)+1

    def genrscript(self):
        with open(self.outname, 'r') as inf, open('out.r','w') as outf:
            outf.write("library(limma)\n"
                      "library(affy)\n")
            outf.write("mydata<-read.table('%s',header = T, row.names = 1, fill = T)\n"%self.outname)
            outf.write("mymatrix <- as.matrix(mydata)\n")
            outf.write("design <- model.matrix(~ 0+factor(c(%s)))\n"%str(self.fill).strip("[]"))

            colnames = ["group"]*(int(self.group)+1)
            for i in range(int(self.group)+1):
                colnames[i] = colnames[i]+str(i+1)
            outf.write("colnames(design) <- c(%s)\n"%str(colnames).strip('[]'))
            outf.write("fit <- lmFit(mymatrix, design)\n")

            contrast = ["group"+str(int(self.group)+1)+"-"+"group"]*int(self.group)
            for i in range(int(self.group)):
                contrast[i] = (contrast[i]+str(i+1))
            outf.write("contrast.matrix <- makeContrasts(%s, levels=design)\n"%','.join(contrast).strip("[]"))
            outf.write('fit2 <- contrasts.fit(fit, contrast.matrix)\n'
                       'fit2 <- eBayes(fit2)\n'
                       'results <- decideTests(fit2)\n'
                       'write.fit(fit2, results, "outputdata.txt", adjust="BH")\n')






# Testing.
project= Pipeline("GSE23006_series_matrix.txt")
project.parser()
project.datainfo()
project.fillline()
project.genrscript()