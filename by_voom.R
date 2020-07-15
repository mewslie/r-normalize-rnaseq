#manually enter paths, files, file attributes etc
file.index <- file.path("C:\\File\\Path\\Here", "filename.csv") #sequence counts
experiment <- data.frame(sample=c("BF1","BF2","AF1","AF2","BM1","BM2","AM1","AM2"), group=rep(rep(c("Before","After"),each=2),2), gender=rep(c("Female","Male"), each=4))
experiment$group <- relevel(as.factor(experiment$group), ref="Before")
experiment$gender <- relevel(as.factor(experiment$gender), ref="Female")
design <- model.matrix(~ experiment$group * experiment$gender) #last factor entered in the formula should be the main factor of interest
colnames(design) <- c("Int","BvsA","MvsF","GroupVGender")
design.group.size <- 2 #no. of samples in comparison group eg GroupVGender
low.count.threshold <- 1 #no. of counts per million that a gene is required to have per sample to stay in this analysis

#load input data
counts.dat <- read.csv(file.index, row.names=1)

#remove sequences with too many low counts
keep <- rowSums(cpm(counts.dat) > low.count.threshold)>=design.group.size; table(keep)
keep.counts <- counts.all[keep,]
###graph counts.keep###
keep.counts.lcpm <- cpm(keep.counts, log=TRUE)
boxplot(keep.counts.lcpm, xlab="", ylab="Log2 counts per million (CPM)", las=2, outline=FALSE, main="logCPM (unnormalised counts)")
abline(h=median(keep.counts.lcpm),col="blue") #mark median logCPM
