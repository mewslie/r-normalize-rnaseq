#manually enter paths, files, file attributes etc
file.index <- file.path("C:\\File\\Path\\Here", "filename.csv") #sequence counts
experiment <- data.frame(sample = c("BF1","BF2","AF1","AF2","BM1","BM2","AM1","AM2"), group = rep(rep(c("Before","After"), each = 2), 2), gender = rep(c("Female","Male"), each = 4))
experiment$group <- relevel(as.factor(experiment$group), ref = "Before")
experiment$gender <- relevel(as.factor(experiment$gender), ref = "Female")
design.group.size <- 2 #no. of samples in comparison group eg GroupVGender
low.count.threshold <- 1 #no. of counts per million that a gene is required to have per sample to stay in this analysis

#load input data
counts.dat <- read.csv(file.index, row.names = 1)

#remove sequences with too many low counts
library("edgeR")
keep <- rowSums(cpm(counts.dat) > low.count.threshold) >= design.group.size; table(keep)
counts.keep <- counts.dat[keep,]
###graph counts.keep###
counts.keep.lcpm <- cpm(counts.keep, log = TRUE)
boxplot(counts.keep.lcpm, xlab = "", ylab = "Log2 counts per million (CPM)", las = 2, outline = FALSE, main = "logCPM (unnormalised counts)")
abline(h = median(counts.keep.lcpm), col = "blue") #mark median logCPM

#use calcNormFactors and voom on counts.keep
counts.dge <- DGEList(counts.keep)
counts.dge <- calcNormFactors(counts.dge, method = "TMM")
#if there is another count table to normalize on eg spike-in sequences, manually change samples$norm.factors
counts.dge$samples$norm.factors <- calcNormFactors("other.counts.table.here", lib.size = counts.dge$samples$lib.size, method = "TMM")
counts.norm <- voom(counts.dge, design, plot = FALSE, normalize.method = "none")
###graph counts on MDS plot###
plotMDS(counts.dge, col = as.numeric(experiment$group))
###graph before and after norm factor calculation###
par(mfrow=c(1,2))
plot.col <- 10 #column of count data to plot
plotMD(cpm(counts.dge, log = TRUE, prior.count = 0.5, normalized.lib.sizes = FALSE), column = plot.col, main = "Unnormalized")
abline(h = 0, col = "grey")
plotMD(counts.norm$E, column = plot.col, main = "Normalized")
abline(h = 0, col = "grey")

#write the counts to file
write.csv(cpm(counts.dge, log = T, prior.count = 0.5), "logcpm_not_normalized.csv")
write.csv(counts.norm$E, "logcpm_normalized.csv")
