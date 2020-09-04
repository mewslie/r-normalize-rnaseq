#manually enter paths, files, file attributes etc
file.index <- file.path("C:\\File\\Path\\Here", "filename.csv") #sequence counts
experiment <- data.frame(sample = c("BF1","BF2","AF1","AF2","BM1","BM2","AM1","AM2"), group = rep(rep(c("Before","After"), each = 2), 2), gender = rep(c("Female","Male"), each = 4))
experiment$group <- relevel(as.factor(experiment$group), ref = "Before")
experiment$gender <- relevel(as.factor(experiment$gender), ref = "Female")
design <- model.matrix(~ experiment$group * experiment$gender) #last factor entered in the formula should be the main factor of interest
colnames(design) <- c("Int","BvsA","MvsF","GroupVGender")
low.count.threshold <- 1 #no. of counts per million that a gene is required to have per sample to stay in this analysis
###graph design###
library("rafalib")
imagemat(design, main = "Model matrix for linear model with interactions", ylab = "Samples", xlab = "Contrasts", xaxt = "n", yaxt = "n")
axis(2, at = 1:nrow(experiment), labels = experiment$SampleName)

#load input data
counts.dat <- read.csv(file.index, row.names = 1)

#remove sequences with too many low counts
library("edgeR")
counts.dge <- DGEList(counts.dat, group = paste(experiment$group, experiment$gender, sep = "."))
counts.keep <- filterByExpr(counts.dge, design, min.count = low.count.threshold)
counts.dge <- counts.dge[counts.keep, , keep.lib.sizes = FALSE]
###graph counts.dge###
counts.lcpm <- cpm(counts.dge, log = TRUE)
boxplot(counts.lcpm, xlab = "", ylab = "Log2 counts per million (CPM)", las = 2, outline = FALSE, main = "logCPM (unnormalised counts)")
abline(h = median(counts.lcpm), col = "blue") #mark median logCPM

#use calcNormFactors and estimateDisp on counts.dge
counts.dge <- calcNormFactors(counts.dge, method = "TMM")
#if there is another count table to normalize on eg spike-in sequences, manually change samples$norm.factors
counts.dge$samples$norm.factors <- calcNormFactors("other.counts.table.here", lib.size = counts.dge$samples$lib.size, method = "TMM")
counts.dge <- estimateDisp(counts.dge, design = design, robust = T)
###graph counts on MDS plot###
plotMDS(counts.dge, col = as.numeric(experiment$group))
###graph before and after norm factor calculation###
par(mfrow=c(1,2))
plot.col <- 10 #column of count data to plot
plotMD(cpm(counts.dge, log = TRUE, prior.count = 0.5, normalized.lib.sizes = FALSE), column = plot.col, main = "Unnormalized")
abline(h = 0, col = "grey")
plotMD(counts.norm$E, column = plot.col, main = "Normalized")
abline(h = 0, col = "grey")
###graph Biological Coefficient of Variation (BVS) plot###
plotBCV(counts.dge)

#write the counts to file
write.csv(cpm(counts.dge, log = T, prior.count = 1), "logcpm_not_normalized.csv")
write.csv(counts.norm$E, "logcpm_normalized.csv")
