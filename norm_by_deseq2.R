#manually enter paths, files, file attributes etc
file.index <- file.path("C:\\File\\Path\\Here", "filename.csv") #sequence counts
experiment <- data.frame(sample = c("BF1","BF2","AF1","AF2","BM1","BM2","AM1","AM2"), group = rep(rep(c("Before","After"), each = 2), 2), gender = rep(c("Female","Male"), each = 4))
design.group.size <- 2 #no. of samples in comparison group eg GroupVGender
low.count.threshold <- 1 #no. of counts per million that a gene is required to have per sample to stay in this analysis

#load input data
counts.dat <- read.csv(file.index, row.names = 1)

#make DESeq2 object
library("DESeq2")
counts.des <- DESeqDataSetFromMatrix(counts.dat, experiment, ~ group + gender + group:gender)
counts.des$Group <- relevel(counts.des$group, "Before") #make sure baseline is control/before
counts.des$Prep <- relevel(counts.des$gender, "Female") #make sure baseline is Female

#remove sequences with too many low counts
keep <- rowSums(counts(counts.des) > low.count.threshold) >= design.group.size; table(keep)
counts.des <- counts.des[keep,]
###graph counts.keep###
counts.log <- normTransform(counts.des)
boxplot(assay(counts.log), xlab="", ylab="Log2(count+1)", las=2, outline=FALSE, main="logC\n(unnormalised counts)")
abline(h=median(assay(counts.log)), col="blue") #mark median log

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
