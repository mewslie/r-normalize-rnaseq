# manually enter paths, files, file attributes etc
rdat.index <- file.path("C:\\File\\Path\\Here", "normalization.Rdata") #Rdata file from previous normalization step; contains raw counts, DGEList object, design matrix
annots.index <- file.path("C:\\File\\Path\\Here", "annotation_table.csv") #addtional annotations for data; one of the columns contains IDs that match the DGEList object
lfc.thres <- 1 #minimum. lfc to accept for DEG
p.thres <- 0.01 #p-value threshold
p.adj <- "BH" #multiple testing adjustment method

# load input data
library("DESeq2")
load(rdat.index)
annot.dat <- read.csv(annots.index)

# run nbinomWaldTest on counts.des; on the last column of model.matrix
counts.des <- nbinomWaldTest(counts.des)

# get top genes
counts.res <- results(counts.des, pAdjustMethod = p.adj, alpha = p.thres))
# use resultsNames to test another contrast
# print summary count of result table
summary(counts.res)
# impose lfc threshold on DEG list
counts.res.lfc <- results(counts.des, pAdjustMethod = "fdr", alpha = p.thres, lfcThreshold = lfc.thres)

# attach annotations to DESeqDataSet (optional)
annots.order <- match(rownames(counts.dge$coefficients), annots.dat$gene.id)
annots.df <- annots.dat[annots.order,]
mcols(counts.des) <- DataFrame(counts.des, annots.df)

# view diagnostic plots for counts.des
###graph MA plot###
library("apeglm")
#shrink log fold change for clearer plot
plot.shrink <- lfcShrink(counts.des, type = "apeglm", coef = 1, res = counts.res)
plotMA(plot.shrink)
###graph dispersion estimates###
plotDispEsts(counts.des)
###graph spread of p-values###
hist(counts.res$pvalue, breaks=20, col="grey50", border="white")
###graph counts of one gene for all samples###
plot.gene <- "gene1"
plotCounts(counts.des, gene = plot.gene, intgroup = c("group","gender"), normalized = TRUE)
###graph volcano plot###
library("EnhancedVolcano")
EnhancedVolcano(counts.res, lab = rownames(counts.res), x = "log2FoldChange", y = "padj")
###use Glimma to do MD plot###
library("Glimma")
glMDPlot(counts.res, status = as.numeric(counts.res$padj < p.thres), counts = counts(counts.ds), groups = targets$group,
  transform = FALSE, samples = colnames(counts.dat), anno = annots.df)
###graph heatmap of top 30 genes###
library("gplots")
get_dist <- function(x,y) {
  topG <- y[order(y$padj),]
  topG <- topG[1:30,]
  expr <- subset(x, row.names(x) %in% row.names(topG))
  expr <- assay(normTransform(expr))
  expr <- t(scale(t(expr)))
  return(expr)
}
plot.dist <- get_dist(counts.des, counts.res, SIMPLIFY = F)
plot.heatmap.col <- colorpanel(100, "blue", "white", "red")
heatmap.2(plot.dist, col = plot.heatmap.col, scale = "none", trace = "none", density.info = "none", dendrogram = "both",
          margin = c(5,9), key = TRUE, key.xlab = "logCPM (scaled)")
###graph PCA of normalized data###
plotPCA(normTransform(counts.des), intgroup=c("group", "gender"))
###graph proportion of low counts vs p-value###
plot.pvals = counts.res$baseMean > 0
plot.quantiles = quantile(counts.res$baseMean[plot.pvals], 0:10/10)
plot.bins = cut(counts.res$baseMean[plot.pvals], qs) #set up the bins
levels(plot.bins) = paste0("~", round(0.5 * qs[-1] + 0.5 * qs[-length(qs)])) #rename the levels of the bins
plot.ratios = tapply(counts.res$padj[plot.pvals], bins, function(p) mean(p < 0.01, na.rm = TRUE)) #calculate the ratio of pvalues > 0.01 for each bin
barplot(ratios, xlab = "Mean normalized count", ylab = "Ratio of small p-values")

# write the result tables to file
counts.out <- counts.res[order(counts.res$padj),])
counts.out <- lapply(counts.out, function(x) {
  index <- match(row.names(x), row.names(mcols(counts.des)))
  x <- cbind(x, mcols(counts.des)[index,])
})
write.csv(counts.out, "topTable.csv")
