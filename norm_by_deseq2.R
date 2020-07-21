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
counts.log <- log2(counts(counts.des, normalized = FALSE)+1)
boxplot(assay(counts.log), xlab = "", ylab = "Log2(count+1)", las = 2, outline = FALSE, main = "logC\n(unnormalised counts)")
abline(h = median(assay(counts.log)), col = "blue") #mark median log

#use estimateSizeFactors on counts.des
counts.des <- estimateSizeFactors(counts.des)
#if there is another count table to normalize on eg spike-in sequences, manually change sizeFactors(counts.des)
sizeFactors(counts.des) <- sizeFactors("other.deseq2.object.here")
###graph before and after norm factor calculation###
counts.norm <- normTransform(counts.des) #default normalized=T, log2
counts.vsd <- vst(counts.des, blind=FALSE) #variance stabilizing transformation
counts.rld <- rlog(counts.des, blind=FALSE) #regularized logarithm
par(mfrow = c(1,4))
boxplot(counts.log, main = "Log2 unnormalized")
abline(h = median(counts.log), col = "red")
boxplot(assay(counts.norm), main = "Log2 Default normalized")
abline(h = median(assay(counts.norm)), col = "red")
boxplot(assay(counts.vsd), main = "VST Default normalized")
abline(h = median(assay(counts.vsd)), col = "red")
boxplot(assay(counts.rld), main = "Rlog Default normalized")
abline(h = median(assay(counts.rld)), col = "red")
###graph PCA of samples###
plotPCA(counts.norm, intgroup=c("group", "gender"))

#write the counts to file
write.csv(counts.log, "log2_not_normalized.csv")
write.csv(assay(counts.norm), "log2_normalized.csv")
