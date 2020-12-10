suppressMessages(library("DESeq2"))
suppressMessages(library("pasilla"))
library(ggplot2)
library(ggfortify)
library(EnhancedVolcano)
library("org.Dm.eg.db")
library(data.table)
library(dplyr)
library(annotate)
source("./geneID_converter.r", local=TRUE)


# Set database to work with
gene_ID_database <- toTable(org.Dm.egFLYBASE)
gene_ID_database_name <- "flybase"

# # Create cts object with non-corrected and unnormalized data 
# cts_df_un <- read.csv("FB_ByCondition_countTable.csv", row.names=1, stringsAsFactors=FALSE)
# cts_mat_un <- as.matrix(cts_df_un, stringsAsFactors = FALSE)
# cols_convert2 <- colnames(cts_df_un)
# for (i in cols_convert2)
# {
#   cts_df_un[[i]] <- as.numeric(cts_df_un[[i]])
# }
# # regex into a list and pass that in to get this
# cts_filtered_adults <- cts_df_un[,c("C8fA.1","C8fA.2","C8fA.3","C8fA.4","C8mA.1","C8mA.2","C8mA.3","C8mA.4",
#                                 "C9fA.1","C9fA.2","C9fA.3","C9fA.4","C9mA.1","C9mA.2","C9mA.3","C9mA.4",
#                                 "C10fA.1","C10fA.2","C10fA.3","C10fA.4","C10mA.1","C10mA.2","C10mA.3","C10mA.4")]

# # first remember the names
# gene_n2 <- rownames(cts_df_un)

# # transpose
# cts_df_un <- as.data.frame(t(cts_df_un))
# colnames(cts_df_un) <- gene_n2

# pca_df_un <- cts_df_un
# pca_df_un <- pca_df_un[ , which(apply(pca_df_un, 2, var) != 0)]
# pca_res_un <- prcomp(pca_df_un, scale. = TRUE)
# png(file="pca_plot_unnorm.png", width=600, height=400)
# autoplot(pca_res_un) + labs(title = "PCA Plot for Un-Normalized and Non-Batched Corrected Read Counts")
# dev.off()
# png(file="pca_plot_unnorm_labeled.png", width=900, height=600, res=125)
# autoplot(pca_res_un, data = cts_condname, colour = 'condition', label = TRUE, label.size = 3) + labs(title = "PCA Plot for Un-Normalized and Non-Batched Corrected Read Counts")
# dev.off()

# # first remember the names
# gene_n2 <- rownames(cts_filtered_adults)

# # transpose
# cts_filtered_adults <- as.data.frame(t(cts_filtered_adults))
# colnames(cts_filtered_adults) <- gene_n2

# # set condition column
# conds_adults = cts_filtered_adults
# conds_adults$condition = rownames(cts_filtered_adults)
# conds_adults$condition <- gsub("\\.\\d", "", conds_adults$condition)
# conds_adults$condition <- gsub("f", "", conds_adults$condition)
# conds_adults$condition <- gsub("m", "", conds_adults$condition)
# conds_adults$condition

# conds_adults$condition[0:2]
# cts_filtered_adults[0:2]
# cts_filtered_adults <- cts_filtered_adults[ , which(apply(cts_filtered_adults, 2, var) != 0)] # find error that caused this...
# pca_filtered_adults <- prcomp(cts_filtered_adults, scale. = TRUE)
# png(file="pca_plot_unnorm_adults.png", width=600, height=400)
# autoplot(pca_filtered_adults) + labs(title = "PCA Plot for Adult Un-Normalized and Non-Batched Corrected Read Counts")
# dev.off()
# png(file="pca_plot_unnorm_labeled_adults.png", width=900, height=600, res=125)
# autoplot(pca_filtered_adults, data = conds_adults, colour = 'condition', label = TRUE, label.size = 3) + labs(title = "PCA Plot for Adult Un-Normalized and Non-Batched Corrected Read Counts")
# dev.off()

# all_samples <- read.csv("metadata_allSamples.csv", row.names=1)
# # url for michael loves answer for pca and normalized counts
# #https://support.bioconductor.org/p/66067/
# #https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
# ddsx <- DESeqDataSetFromMatrix(countData = cts_mat_un, colData = all_samples, design = ~ condition + batch) # + batch
# dds_temp1 <- estimateSizeFactors(ddsx)
# dds_temp2 <- counts(dds_temp1, normalized=TRUE)

# write.csv(dds_temp2, file="normalized_counts.csv")


# cts_filtered_adults <- dds_temp2[,c("C8fA.1","C8fA.2","C8fA.3","C8fA.4","C8mA.1","C8mA.2","C8mA.3","C8mA.4",
#                                 "C9fA.1","C9fA.2","C9fA.3","C9fA.4","C9mA.1","C9mA.2","C9mA.3","C9mA.4",
#                                 "C10fA.1","C10fA.2","C10fA.3","C10fA.4","C10mA.1","C10mA.2","C10mA.3","C10mA.4")]

# # first remember the names
# gene_n22 <- rownames(dds_temp2)

# # transpose
# dds_temp2 <- as.data.frame(t(dds_temp2))
# colnames(dds_temp2) <- gene_n22

# dds_temp2 <- dds_temp2[ , which(apply(dds_temp2, 2, var) != 0)]
# pca_res_deseqnorm <- prcomp(dds_temp2, scale. = TRUE)
# png(file="pca_plot_deseqnorm.png", width=600, height=400)
# autoplot(pca_res_deseqnorm) + labs(title = "PCA Plot for Normalized Read Counts using DESeq2 Method")
# dev.off()
# png(file="pca_plot_deseqnorm_labeled.png", width=900, height=600, res=125)
# autoplot(pca_res_deseqnorm, data = cts_condname, colour = 'condition', label = TRUE, label.size = 3) + labs(title = "PCA Plot for Normalized Read Counts using DESeq2 Method")
# dev.off()

# # first remember the names
# gene_n22 <- rownames(cts_filtered_adults)

# # transpose
# cts_filtered_adults <- as.data.frame(t(cts_filtered_adults))
# colnames(cts_filtered_adults) <- gene_n22

# cts_filtered_adults <- cts_filtered_adults[ , which(apply(cts_filtered_adults, 2, var) != 0)]
# pca_filtered_adults <- prcomp(cts_filtered_adults, scale. = TRUE)
# png(file="pca_plot_deseqnorm_adults.png", width=600, height=400)
# autoplot(pca_filtered_adults) + labs(title = "PCA Plot for Adult Normalized Read Counts using DESeq2 Method")
# dev.off()
# png(file="pca_plot_deseqnorm_labeled_adults.png", width=900, height=600, res=125)
# autoplot(pca_filtered_adults, data = conds_adults, colour = 'condition', label = TRUE, label.size = 3) + labs(title = "PCA Plot for Adult Normalized Read Counts using DESeq2 Method")
# dev.off()


cts <- as.matrix(read.csv("FB_ByCondition_countTable.csv", row.names=1, stringsAsFactors=FALSE), stringsAsFactors=FALSE)
# coldata <- read.csv("metadata.csv")
coldata <- read.csv("metadata_allSamples.csv", row.names=1)
sample_nms <- rownames(coldata)
rows_select <- grep("^(C4e|C6e)",sample_nms,value=TRUE)
coldata <- coldata[c(rows_select),]
coldata <- coldata[order(coldata$ID),]
rownames(coldata) <- coldata$ID
coldata$ID <- NULL
head(cts,2)
head(coldata,2)
nms <- colnames(cts)
cols_select <- grep("^(C4e|C6e)",nms,value=TRUE)
cts_filtered <- cts[,c(cols_select)]
for (i in cols_select)
{
  cts_filtered[,i] <- as.integer(cts_filtered[,i])
}
# cts_filtered <- cts[,c("C4e.1","C4e.2","C4e.3","C4e.4","C6e.1","C6e.2","C6e.3","C6e.4")]
# cts_filtered <- transform(cts_filtered, C4e.1 = as.integer(C4e.1))
# cts_filtered <- transform(cts_filtered, C4e.2 = as.integer(C4e.2))
# cts_filtered <- transform(cts_filtered, C4e.3 = as.integer(C4e.3))
# cts_filtered <- transform(cts_filtered, C4e.4 = as.integer(C4e.4))
# cts_filtered <- transform(cts_filtered, C6e.1 = as.integer(C6e.1))
# cts_filtered <- transform(cts_filtered, C6e.2 = as.integer(C6e.2))
# cts_filtered <- transform(cts_filtered, C6e.3 = as.integer(C6e.3))
# cts_filtered <- transform(cts_filtered, C6e.4 = as.integer(C6e.4))
head(cts_filtered,2)
all(rownames(coldata) %in% colnames(cts_filtered))
all(rownames(coldata) == colnames(cts_filtered))
cts_filtered <- cts_filtered[, rownames(coldata)]
all(rownames(coldata) == colnames(cts_filtered))
ncol(cts_filtered)
colnames(cts_filtered)
nrow(coldata)
rownames(coldata)
#sapply(cts_filtered, class)
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
resNorm_no_padj <- lfcShrink(dds, coef=1, type="normal") # coef=1 is elav_gfp vs rnai 
resNorm <- resNorm_no_padj[which(resNorm_no_padj$padj < 0.01),]
write.csv(res, file = 'deseq2_results_pvalandnorm1.csv')


# pdf("volcano_plot_isaacdeseq2_pvalfiltered_unnorm.pdf", width = 15, height = 15)
# res1 = results(dds, contrast = c('condition', 'elavGRFP_e', 'cRNAi_e'))
# resNorm_no_padj1 <- lfcShrink(dds, contrast = c('condition', 'elavGRFP_e', 'cRNAi_e'), res=res1, type="normal") # coef=1 is elav_gfp vs rnai 
# resNorm1 <- resNorm_no_padj1[which(resNorm_no_padj1$padj < 0.01),]
# # Adding gene symbol and placing it in the front for no and normal shrinkage matricies
# ids.type  <- gene_ID_database_name
# idsN <- rownames(resNorm1)
# resNorm1['gene_id'] <- rownames(resNorm1)
# resNorm1$gene_name <- as.vector(get.symbolIDsDm(idsN,ids.type))
# res_sym_front <- as.data.frame(resNorm1) %>% dplyr::select(gene_name, gene_id, everything())  
# EnhancedVolcano(res_sym_front, 
#                 lab = res_sym_front$gene_name, 
#                 x = 'log2FoldChange', 
#                 y = 'padj', 
#                 xlim = c(-5, 5),
#                 ylim = c(0, 35),
#                 pCutoff = 0.01,
#                 FCcutoff = 0,
#                 pointSize = 3.0,
#                 labSize = 6.0,
#                 title = "Volcano Plot for Un-Normalized Read Counts Adjusted by DESeq2 Normalization")
# dev.off()
