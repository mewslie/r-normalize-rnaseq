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

# run glmQLFit on counts.dge
counts.dge <- glmQLFit(counts.dge, design)
# run quasi-likelihood F-test on counts.dge; on the last column of model.matrix or your custom contrast (check model colnames for contrast names) with makeContrasts
counts.glm <- glmQLFTest(counts.dge)
# impose lfc threshold on DEG list
counts.glm.lfc <- glmTreat(counts.dge, lfc = lfc.thres)
# get top genes
counts.glm.top <- topTags(counts.glm, n = Inf, sort.by = "p", adjust.method = p.adj)
# print summary count of test results
summary(decideTests(counts.glm, adjust.method = p.adj, p.value = p.thres))
# attach annotations to result table
annots.order <- match(rownames(counts.dge$coefficients), annots.dat$gene.id)
annots.df <- annots.dat[annots.order,]
counts.glm$genes <- annots.df

# view diagnostic plots for counts.glm
###graph MD plot###
plotMD(counts.glm, hl.col = c("red","blue"))
###graph smear plot aka edgeR-version of mean-difference plot###
plot.tags <- decideTestsDGE(counts.glm, adjust.method = p.adj, p.value = p.thres
plotSmear(counts.glm, de.tags = rownames(plot.tags))
###graph quasi likelihood dispersion plot###
plotQLDisp(counts.glm)
###plot heatmap###
library("gplots")
heatmap.genes <- 30 #plot top 30 genes
lfc.norm <- cpm(counts.dge, log = T, prior.count = 0.5, normalized.lib.sizes = T)
heatmap.df <- merge(lfc.norm, counts.glm.top[1:heatmap.genes, c(1:2,9,10)], by = 0)
heatmap.df <- t(scale(t(heatmap.df[,2:13])))
heatmap.col <- colorpanel(100, "blue", "white", "red")
heatmap.2(heatmap.df, col = heatmap.col, scale = "none", trace = "none", density.info = "none", dendrogram = "both",
          margin = c(5,7), key = TRUE, key.xlab = "logCPM (scaled)")

#write the result tables to file
write.csv(counts.glm.top, "topTable.csv")
