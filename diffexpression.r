suppressMessages(library("DESeq2"))
suppressMessages(library("pasilla"))
# library(data.table)
library(ggplot2)
library(ggfortify)
# cts <- as.matrix(read.csv("afterbatchtable_geneNameAdded.csv", row.names=1, stringsAsFactors=FALSE), stringsAsFactors=FALSE)
cts_df <- read.csv("afterbatchtable_geneNameAdded.csv", row.names=1, stringsAsFactors=FALSE)
cts_mat <- as.matrix(cts_df, stringsAsFactors = FALSE)
all_samples <- read.csv("metadata_allSamples.csv", row.names=1)
# metaData$time <- as.factor(metaData$time)
head(cts_df,2)
cts_df <- cts_df[1:68]
cols_convert <- colnames(cts_df)
print(cols_convert)
for (i in cols_convert)
{
    #cts_df <- transform(cts_df, i = as.integer(i))
#   cts_df$i <- as.numeric(cts_df$i)
  cts_df[[i]] <- as.numeric(cts_df[[i]])

}
# cts_df$i <- NULL
sapply(cts_df, class)

# first remember the names
gene_n <- rownames(cts_df)

# transpose
cts_df <- as.data.frame(t(cts_df))
colnames(cts_df) <- gene_n
# df.aree$myfactor <- factor(row.names(df.aree))

# str(df.aree) # Check the column types
cts_df[c(1:5),c(1:5)]
cts_condname = cts_df
cts_condname$condition = rownames(cts_df)
cts_condname$condition <- gsub("\\.\\d", "", cts_condname$condition)
cts_condname$condition <- gsub("f", "", cts_condname$condition)
cts_condname$condition <- gsub("m", "", cts_condname$condition)
cts_condname$condition
# cts_condname$condition[cts_condname$condition == "C3e.1"] <- "b"
# cts_condname$condition[cts_condname$condition == "B"] <- "b"
# cts_condname$condition[cts_condname$condition == "B"] <- "b"
# cts_condname$condition[cts_condname$condition == "B"] <- "b"
# cts_condname$condition[cts_condname$condition == "B"] <- "b"
# cts_condname$condition[cts_condname$condition == "B"] <- "b"
# cts_condname$condition[cts_condname$condition == "B"] <- "b"
# cts_condname$condition[cts_condname$condition == "B"] <- "b"
# cts_condname$condition[cts_condname$condition == "B"] <- "b"
# cts_condname$condition[cts_condname$condition == "B"] <- "b"
# cts_condname$condition[cts_condname$condition == "B"] <- "b"
# cts_condname$condition[cts_condname$condition == "B"] <- "b"
# cts_condname$condition[cts_condname$condition == "B"] <- "b"

# # head(cts_df,2)
# sum(is.na(cts_df))
# pca_df <- cts_df
# pca_res <- prcomp(pca_df, scale. = TRUE)
# png(file="pca_plot.png")
# autoplot(pca_res)
# dev.off()
# png(file="pca_plot_labeled.png", width=900, height=600, res=150)
# autoplot(pca_res, data = cts_condname, colour = 'condition', label = TRUE, label.size = 3)
# dev.off()

#TODO: create cts object with non-corrected and normalized data 
cts_df_un <- read.csv("FB_ByCondition_countTable.csv", row.names=1, stringsAsFactors=FALSE)
cts_mat_un <- as.matrix(cts_df_un, stringsAsFactors = FALSE)
head(cts_df_un,2)
cols_convert2 <- colnames(cts_df_un)
print(cols_convert2)
for (i in cols_convert2)
{
    #cts_df <- transform(cts_df, i = as.integer(i))
#   cts_df$i <- as.numeric(cts_df$i)
  cts_df_un[[i]] <- as.numeric(cts_df_un[[i]])

}
# first remember the names
gene_n2 <- rownames(cts_df_un)

# transpose
cts_df_un <- as.data.frame(t(cts_df_un))
colnames(cts_df_un) <- gene_n2
# print(gene_n2)
# df.aree$myfactor <- factor(row.names(df.aree))

# str(df.aree) # Check the column types
cts_df_un[c(1:5),c(1:5)]
# cts_condname = cts_df
# cts_condname$condition = rownames(cts_df)
# cts_condname$condition <- gsub("\\.\\d", "", cts_condname$condition)
# cts_condname$condition <- gsub("f", "", cts_condname$condition)
# cts_condname$condition <- gsub("m", "", cts_condname$condition)
# cts_condname$condition
# pca_df_un <- cts_df_un
# pca_df_un <- pca_df_un[ , which(apply(pca_df_un, 2, var) != 0)]
# pca_res_un <- prcomp(pca_df_un, scale. = TRUE)
# png(file="pca_plot_unnorm.png", width=600, height=400)
# autoplot(pca_res_un)
# dev.off()
# png(file="pca_plot_unnorm_labeled.png", width=900, height=600, res=150)
# autoplot(pca_res_un, data = cts_condname, colour = 'condition', label = TRUE, label.size = 3)
# dev.off()

cts <- as.matrix(read.csv("afterbatchtable_geneNameAdded.csv", row.names=1, stringsAsFactors=FALSE), stringsAsFactors=FALSE)
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
# dds_temp1 <- estimateSizeFactors(dds)
# dds_temp2 <- counts(dds_temp1, normalized=TRUE)
# dds_temp2
dds <- DESeq(dds, test="LRT", reduced = ~1) #batch
res <- results(dds)
res
resultsNames(dds)
# write.csv(res, file = 'deseq2_results1.csv')
# Filtering based on adjusted p-value (i.e. padj)
res_no_padj <- results(dds)
res <- res_no_padj[which(res_no_padj$padj < 0.01),]
# write.csv(res, file = 'deseq2_results_pvalfiltered1.csv')

# Performing normal shrinkage transformation

#TODO: AC figure out best shrinkage
resNorm_no_padj <- lfcShrink(dds, coef=1, type="normal") # coef=1 is elav_gfp vs rnai 
resNorm <- resNorm_no_padj[which(resNorm_no_padj$padj < 0.01),]
# write.csv(res, file = 'deseq2_results_pvalandnorm1.csv')



# cts_df_un <- as.matrix(read.csv("FB_ByCondition_countTable.csv", row.names=1, stringsAsFactors=FALSE), stringsAsFactors=FALSE)
# # head(cts_df_un, 2)
# # xasx <- sapply(cts_df_un, class)
# # xasx[0:2]
# # cols_convert2 <- colnames(cts_df_un)
# # for (i in cols_convert2)
# # {
# #   cts_df_un[[i]] <- as.integer(cts_df_un[[i]])

# # }
# ddsx <- DESeqDataSetFromMatrix(countData = cts_df_un, colData = all_samples, design = ~ condition) # + batch
# dds_temp1 <- estimateSizeFactors(ddsx)
# dds_temp2 <- counts(dds_temp1, normalized=TRUE)
# # dds_temp2[0:2]
# write.csv(dds_temp2, file="normalized_counts.csv")
# # first remember the names
# gene_n22 <- rownames(dds_temp2)

# # transpose
# dds_temp2 <- as.data.frame(t(dds_temp2))
# colnames(dds_temp2) <- gene_n22
# # print(gene_n2)
# # df.aree$myfactor <- factor(row.names(df.aree))

# # str(df.aree) # Check the column types
# dds_temp2[c(1:5),c(1:5)]
# # cts_condname = cts_df
# # cts_condname$condition = rownames(cts_df)
# # cts_condname$condition <- gsub("\\.\\d", "", cts_condname$condition)
# # cts_condname$condition <- gsub("f", "", cts_condname$condition)
# # cts_condname$condition <- gsub("m", "", cts_condname$condition)
# # cts_condname$condition
# dds_temp2 <- dds_temp2[ , which(apply(dds_temp2, 2, var) != 0)]
# pca_res_deseqnorm <- prcomp(dds_temp2, scale. = TRUE)
# png(file="pca_plot_deseqnorm.png", width=600, height=400)
# autoplot(pca_res_deseqnorm)
# dev.off()
# png(file="pca_plot_deseqnorm_labeled.png", width=900, height=600, res=150)
# autoplot(pca_res_deseqnorm, data = cts_condname, colour = 'condition', label = TRUE, label.size = 3)
# dev.off()

library(EnhancedVolcano)

# # png(file="volcano_plot_deseq2v3.png", width=1000, height=1000, res=200)
# pdf("volcano_plot_deseq2_new.pdf", width = 15, height = 15)
# res1 = results(dds, contrast = c('condition', 'elavGRFP_e', 'cRNAi_e'))
# resNorm_no_padj1 <- lfcShrink(dds, contrast = c('condition', 'elavGRFP_e', 'cRNAi_e'), res=res1, type="normal") # coef=1 is elav_gfp vs rnai 
# EnhancedVolcano(resNorm_no_padj1, 
#                 lab = rownames(resNorm_no_padj1), 
#                 x = 'log2FoldChange', 
#                 y = 'padj', 
#                 xlim = c(-5, 5),
#                 ylim = c(0, 35),
#                 pCutoff = 0.01,
#                 FCcutoff = 0.5,
#                 pointSize = 3.0,
#                 labSize = 3.0)
# dev.off()

# pdf("volcano_plot_isaacdeseq2_pvalfiltered.pdf", width = 15, height = 15)
# res1 = results(dds, contrast = c('condition', 'elavGRFP_e', 'cRNAi_e'))
# resNorm_no_padj1 <- lfcShrink(dds, contrast = c('condition', 'elavGRFP_e', 'cRNAi_e'), res=res1, type="normal") # coef=1 is elav_gfp vs rnai 
# resNorm1 <- resNorm_no_padj1[which(resNorm_no_padj1$padj < 0.01),]
# EnhancedVolcano(resNorm_no_padj1, 
#                 lab = rownames(resNorm_no_padj1), 
#                 x = 'log2FoldChange', 
#                 y = 'padj', 
#                 xlim = c(-5, 5),
#                 ylim = c(0, 35),
#                 pCutoff = 0.01,
#                 FCcutoff = 0,
#                 pointSize = 3.0,
#                 labSize = 3.0)
# dev.off()

# pdf("volcano_plot_ashleydebrowser_pvalfiltered.pdf", width = 15, height = 15)
# resNorm_no_padj11 <- as.data.frame(read.csv("combined_clampi_e.csv", row.names=1, stringsAsFactors=FALSE))
# # cols_convert_p <- colnames(resNorm_no_padj11)
# # for (i in cols_convert_p)
# # {
# #   resNorm_no_padj11[[i]] <- as.numeric(resNorm_no_padj11[[i]])
# # }
# EnhancedVolcano(resNorm_no_padj11, 
#                 lab = rownames(resNorm_no_padj11), 
#                 x = 'log2FoldChange', 
#                 y = 'padj',
#                 pCutoff = 0.01,
#                 FCcutoff = 0,
#                 pointSize = 3.0,
#                 labSize = 3.0)
# dev.off()



cts <- as.matrix(read.csv("FB_ByCondition_countTable.csv", row.names=1, stringsAsFactors=FALSE), stringsAsFactors=FALSE)
coldata <- read.csv("metadata.csv")
coldata <- coldata[order(coldata$ID),]
rownames(coldata) <- coldata$ID
coldata$ID <- NULL
head(cts,2)
head(coldata,2)

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
# write.csv(res, file = 'deseq2_results1.csv')
# Filtering based on adjusted p-value (i.e. padj)
res_no_padj <- results(dds)
res <- res_no_padj[which(res_no_padj$padj < 0.01),]
# write.csv(res, file = 'deseq2_results_pvalfiltered1.csv')

# Performing normal shrinkage transformation

#TODO: AC figure out best shrinkage
resNorm_no_padj <- lfcShrink(dds, coef=1, type="normal") # coef=1 is elav_gfp vs rnai 
resNorm <- resNorm_no_padj[which(resNorm_no_padj$padj < 0.01),]
# write.csv(res, file = 'deseq2_results_pvalandnorm1.csv')


pdf("volcano_plot_isaacdeseq2_pvalfiltered_unnorm.pdf", width = 15, height = 15)
res1 = results(dds, contrast = c('condition', 'elavGRFP_e', 'cRNAi_e'))
resNorm_no_padj1 <- lfcShrink(dds, contrast = c('condition', 'elavGRFP_e', 'cRNAi_e'), res=res1, type="normal") # coef=1 is elav_gfp vs rnai 
resNorm1 <- resNorm_no_padj1[which(resNorm_no_padj1$padj < 0.01),]
EnhancedVolcano(resNorm_no_padj1, 
                lab = rownames(resNorm_no_padj1), 
                x = 'log2FoldChange', 
                y = 'padj', 
                xlim = c(-5, 5),
                ylim = c(0, 35),
                pCutoff = 0.01,
                FCcutoff = 0,
                pointSize = 3.0,
                labSize = 3.0)
dev.off()
