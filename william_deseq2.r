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
dds <- DESeqDataSetFromMatrix(countData = cts_filtered, colData = coldata, design = ~ 0 + condition + batch)
dds
dds <- DESeq(dds, test="LRT", reduced = ~batch)
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


# set up rnai vs. control ... over the time course 
design(ddsMF) <- formula(~ time + condition + batch + sex + time:condition)
ddsMF <- DESeq(ddsMF)
resMF <- results(ddsMF)
resMFType <- results(ddsMF, #reference - control is reference
                     contrast=c("sex", "female", "male")) # female (numerator), male (denominator) reference)

e_f_rnai vs. control - 10 DEG
L3_f_rnai vs. control - 20 DEG
A_f_rnai vs. control (*) - 10 DEG
Total DEG = 40 

# look at L3 and A together and use embryo as a reference
dds <- embryo + control

metadata - samples, stage, sex, treatment, batch   
                    embryo  
                    L3
                    A 

# levels(met$treatment) <- factor(met$treatment, levels(c("REF (control)", "rnai")))
# DO FOR EVERY COLUMN
# or a way around it
# ~ 0 +  <- removing reference!!!!!!!
# treatment - elavGFP, rnai
# time - e, L3, A 

# interaction terms in the linear model
# 1)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ 0 + batch + treatment + stage + sex + treatment:time, treatment:stage, treatment:sex, batch:treatment:stage) # examine effect of condition over time
# 2)
dds <- DESeq(dds, test="LRT", reduced= ~ treatment, treatment:time, treatment:time, treatment:stage, treatment:sex, batch:treatment:stage) # whatever is important 
# 3)
res <- resultsNames(dds)
# choose the column names - res <- resultsNames(dds)

# filter the larger results to look for what you want to compare
#?coefficient
#?intercept
# t-test (x dif than y)
# I dont have to have a reference in contrast!
# do ?contrasts
# rerun 1
# 4)
results(dds, contrast = c("treatment:time") #, "rnai:e", "rnai:L3", "rnai:A", "elavGFP:e", "elavGFP:L3", "elavGFP:A"))# denom is reference!

# rerun 1
# another contrast
results(dds, contrast = c("treatment", "rnai", "elavGFP"))# denom is reference!

# columns in results : control:embryo, rnai:embryo

# find the columns of interest and specify that as a comparison
# issue with referecez

#goal: effect over time
