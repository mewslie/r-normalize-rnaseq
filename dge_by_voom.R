# manually enter paths, files, file attributes etc
rdat.index <- file.path("C:\\File\\Path\\Here", "normalization.Rdata") #Rdata file from previous normalization step; contains raw counts, DGEList object, design matrix
annots.index <- file.path("C:\\File\\Path\\Here", "annotation_table.csv") #addtional annotations for data; one of the columns contains IDs that match the DGEList object
lfc.thres <- 1 #minimum. lfc to accept for DEG
p.thres <- 0.01 #p-value threshold
p.adj <- "BH" #multiple testing adjustment method

# load input data
library("edgeR")
load(rdat.index)
annot.dat <- read.csv(annots.index)

# run lmfit
counts.lm <- lmFit(counts.norm, design)
# (optional) make a contrast for a specific comparison
test.contrast <- makeContrasts(groupFemale - groupMale, levels = colnames(coef(countslm)))
test.contrast.fit <- contrasts.fit(counts.lm, test.contrast)
# run eBayes on fitted model
test.fit <- eBayes(counts.lm)
# check number of DEGs
summary(decideTests(test.fit, adjust.method = p.adj, p.value = p.thres)
# get top table
top.table <- topTable(test.fit, sort.by = "P", n = Inf)
# print number of DEGs
length(which(top.table$adj.P.Val < 0.05))
# (optional) impose lfc threshold to filter down DEG
test.fit.lfc <- treat(test.fit, lfc = lfc.thres)

# diagnostic plots
###graph MD plot###
test.res <- decideTests(counts.lm, adjust.method = p.adj, p.value = p.thres)
plotMD(counts.lm, coef = 3, status = test.res[,"femalevsmale"], values = c("1","-1"), hl.col = c("green","red"))
###graph volcano plot###
volcanoplot(counts.lm, coef=3, highlight = 100, names = counts.lm$genes$Name)

# write top results to file
write.csv(top.table,"topTable.csv")
