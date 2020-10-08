suppressMessages(library("DESeq2"))
suppressMessages(library("pasilla"))
cts <- as.matrix(read.csv("afterbatchtable_geneNameAdded.csv", row.names=1, stringsAsFactors=FALSE), stringsAsFactors=FALSE)

#TODO: create cts object with non-corrected and normalized data 

coldata <- read.csv("metadata.csv")
coldata <- coldata[order(coldata$ID),]
rownames(coldata) <- coldata$ID
coldata$ID <- NULL
head(cts,2)
head(coldata,2)

#TODO: PCA plot - label dots

cts_filtered <- cts[,c("C4e.1","C4e.2","C4e.3","C4e.4","C6e.1","C6e.2","C6e.3","C6e.4")]
cts_filtered <- transform(cts_filtered, C4e.1 = as.integer(C4e.1))
cts_filtered <- transform(cts_filtered, C4e.2 = as.integer(C4e.2))
cts_filtered <- transform(cts_filtered, C4e.3 = as.integer(C4e.3))
cts_filtered <- transform(cts_filtered, C4e.4 = as.integer(C4e.4))
cts_filtered <- transform(cts_filtered, C6e.1 = as.integer(C6e.1))
cts_filtered <- transform(cts_filtered, C6e.2 = as.integer(C6e.2))
cts_filtered <- transform(cts_filtered, C6e.3 = as.integer(C6e.3))
cts_filtered <- transform(cts_filtered, C6e.4 = as.integer(C6e.4))
head(cts_filtered,2)
all(rownames(coldata) %in% colnames(cts_filtered))
all(rownames(coldata) == colnames(cts_filtered))
cts_filtered <- cts_filtered[, rownames(coldata)]
all(rownames(coldata) == colnames(cts_filtered))
ncol(cts_filtered)
colnames(cts_filtered)
nrow(coldata)
rownames(coldata)
sapply(cts_filtered, class)
dds <- DESeqDataSetFromMatrix(countData = cts_filtered, colData = coldata, design = ~ condition) # + batch
dds
dds <- DESeq(dds, test="LRT", reduced = ~1) #batch
res <- results(dds)
res
resultsNames(dds)
write.csv(res, file = 'deseq2_results1.csv')
# Filtering based on adjusted p-value (i.e. padj)
res_no_padj <- results(dds)
res <- res_no_padj[which(res_no_padj$padj < 0.01),]
write.csv(res, file = 'deseq2_results_pvalfiltered1.csv')

# Performing normal shrinkage transformation

#TODO: AC figure out best shrinkage
resNorm_no_padj <- lfcShrink(dds, coef=1, type="normal") # coef=1 is elav_gfp vs rnai 
resNorm <- resNorm_no_padj[which(resNorm_no_padj$padj < 0.01),]
write.csv(res, file = 'deseq2_results_pvalandnorm1.csv')



