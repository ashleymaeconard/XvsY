# DESeq2.r
# Purpose: run DESeq2 in case vs. control conditions
# Input: 1) metadata file, 2) count matrix file, 3) desired output directory, 4) condition name, and are there 5) batch effects or 6) time components for the experiment set at 7) adj. p-value threshold  
# Output: 1) MA plots (normal and no shrinkage), 2) venn diagram of LRT vs. Wald test (no shrinkage), 3) clustermap differentially expressed gene .csv (normal and no shrinkage), 4) DESeq2 differentially expressed gene output .csv (normal and no shrinkage)
# Last Mod: Nov. 19, 2010
# Ashley Conard

#Input arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Type:DESeq2.r /FULL/PATH/TO/METDATA_FILE/ /FULL/PATH/TO/COUNT_MATRIX_FILE/ (not normalized and corrected) /FULL/PATH/TO/OUTPUTDIR/ CONDITION (e.g. insulin_stim) BATCH_EFFECT (0 or 1) TIME_COURSE (0 or 1) ADJP_THRESH", call.=FALSE)
} else if (length(args) == 7) {
    cat("Passed in:", args,"\n")
} else{
    stop("Pass in 7 arguments. Type:DESes2.r 1) /FULL/PATH/TO/METDATA_FILE/ 2) /FULL/PATH/TO/COUNT_MATRIX_FILE/ (not normalized and corrected) 3) /FULL/PATH/TO/OUTPUTDIR/ 4) CONDITION (e.g. insulin_stim) 5) BATCH_EFFECT (0 or 1) 6) TIME_COURSE (0 or 1) ADJP_THRESH")
}

METDATA <- args[1] 
COUNTDATA <- args[2]
OUTPUTDIR <- args[3] 
CONDITION <- args[4]
BATCH_EFFECT <- as.numeric(args[5]) # 0 or 1, where 1 sets the reduced model to calculate the effect of time, controlling for the effect of batch (in reduced model)
TC <- as.numeric(args[6]) # 0 or 1, if time course then 1
PVAL_THRESH <- as.numeric(args[7])

print(METDATA)
print(COUNTDATA)
print(CONDITION)

# Libraries
library("DESeq2")# make sure 1.22
library(dplyr)
library(annotate)
library("org.Dm.eg.db")
source("./clampvsmsl2/geneID_converter.r", local=TRUE)
library(data.table)

# Main function
main <- function(){
    
    print("CONDITION")
    print(CONDITION)
    
    # Set database to work with
    gene_ID_database <- toTable(org.Dm.egFLYBASE)
    gene_ID_database_name <- "flybase"

    # Create or set output directory
    CONDIT_DIR <- paste(CONDITION,"results",sep="_")
    FULL_OUTDIR <- paste(OUTPUTDIR,CONDIT_DIR,"deseq2",sep="/")
    if (!dir.exists(CONDIT_DIR)){
        dir.create(CONDIT_DIR)
    } else {
        print("Results directory already exists.")
    }
    if (!dir.exists(FULL_OUTDIR)){
        dir.create(FULL_OUTDIR)
    } else {
        print("deseq2 directory already exists.")
    }
    cat("Output directory: ", FULL_OUTDIR)

    # Load data and metadata
    # Import count data
    countData1 <- read.csv(COUNTDATA, header=TRUE, sep=",")
    col_name = "ID"
    # Import metadata
    metaData <- read.csv(METDATA, header=TRUE, sep=",")
    col_for_index="ID"

    # Format metadata
    keep1 <- as.vector(metaData[[col_for_index]]) # conditions to keep
    rownames(metaData) <- metaData$ID # set rownames to be condition labels
    drop<-c(col_for_index)
    metaData = metaData[,!(names(metaData) %in% drop)]
    metaData$time <- as.factor(metaData$time) # use factor for Time (which is all integers)
    metaData$condition <- as.factor(metaData$condition)
    metaData$batch <- as.factor(metaData$batch)
    metaData <- metaData[ order(row.names(metaData)), ]
    
    # Checking to see if case vs. control or just case or control
    if (dim(table(metaData$condition)) == 1){ 
        CASEvsCONT <- 0  
    } else{
        CASEvsCONT <- 1 # yes this is case vs. control of some sort
        
        # Checking case vs. control (i.e. do all timepoints have a control or is it one control for all other timepoints)
        freq_case_contr = data.frame(table(metaData$condition))
    }
    print(head(metaData))

    # Make sure all Flybase ID names are unique
    n_occur <- data.frame(table(countData1[col_name]))
    non_unique = n_occur[n_occur$Freq > 1,]
    if(length(non_unique) == 0){
        print("ERROR - COL. NAMES (FLYBASE) NOT ALL UNIQUE")
        print(nrow(countData1))
        print(non_unique)
        print(countData1[countData1[col_name] %in% n_occur$Var1[n_occur$Freq > 1],])
        quit()
    }
    rownames(countData1) <- countData1[[col_name]]# use the Flybase name as unique ID
    countData1 <- countData1[, !(names(countData1) %in% col_name)] # remove Flybase name from columns

    # Subsetting merged_htseq based on metadata file
    ids_to_keep <- c(rownames(metaData))
    print("ids_to_keep")
    print(ids_to_keep)
    col.num <- which(colnames(countData1) %in% ids_to_keep)
    countData <- countData1[,col.num]
    countData <- countData[ , order(names(countData))]
    print(head(countData))

    # Check that all row names for metaData match column names for countData
    if(!all(rownames(metaData) %in% colnames(countData))){
        write("ERROR: row names for metaData do not match column names for countData.", stderr())
    }
    
    # BATCH_EFFECT & TC & CASEvsCONT

    # Create DESeq2 object and run DESeq2
    if(BATCH_EFFECT & !TC & CASEvsCONT){ # batch effect and condition
        write("Batch effect and case vs. control", stderr())
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ batch + condition) # examine effect of condition
        dds <- DESeq(dds, test="Wald", betaPrior=TRUE, reduced= ~ batch) ## add Wald test here
    } else if(BATCH_EFFECT & TC & !CASEvsCONT){ # batch effect and timecourse
        write("Batch effect and timecourse", stderr())
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ batch + time)# examine effect of time
        dds <- DESeq(dds, test="LRT", reduced= ~ batch)
    # 1 variable
    } else if(!BATCH_EFFECT & !TC & CASEvsCONT){ # condition (normally what DESeq2 is used for - case vs. control, no time)
        write("Case vs. control", stderr())
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ condition)
        dds <- DESeq(dds, test="Wald", betaPrior=TRUE, reduced= ~1)
    } else if(!BATCH_EFFECT & TC & !CASEvsCONT){ # condition (normally what DESeq2 is used for - case vs. control, no time)
        write("Timecourse", stderr())
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ time)
        dds <- DESeq(dds, test="LRT", reduced= ~1)
    } else if(BATCH_EFFECT & !TC & !CASEvsCONT){ # condition (normally what DESeq2 is used for - case vs. control, no time)
        write("Batch effect", stderr())
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ batch)
        dds <- DESeq(dds, test="LRT", reduced= ~1)
    } else{
        write("ERROR: TIMEOR processes datasets with columns 'time', 'batch', and 'condition'.", stderr())
    }

    # Filtering based on adjusted p-value (i.e. padj)
    res_no_padj <- results(dds)
    res <- res_no_padj[which(res_no_padj$padj < PVAL_THRESH),]
    
    # Performing normal shrinkage transformation
    resNorm_no_padj <- lfcShrink(dds, coef=3, type="normal") # coef=3 is treatment_odor_vs_etoh  
    resNorm <- resNorm_no_padj[which(resNorm_no_padj$padj < PVAL_THRESH),]

    # Saving MA plots before and after shrinkage
    pdf(paste(FULL_OUTDIR,paste('ma_plot_noShrinkage_padj',toString(PVAL_THRESH),'.pdf',sep=""),sep="/"))
    plotMA(res, ylim=c(-4,4), cex=.8)
    abline(h=c(-1,1), col="dodgerblue", lwd=2)
    dev.off()
    pdf(paste(FULL_OUTDIR,paste('ma_plot_normalShrinkage_padj',toString(PVAL_THRESH),'.pdf',sep=""), sep="/"))
    plotMA(resNorm, ylim=c(-4,4), cex=.8)
    abline(h=c(-1,1), col="dodgerblue", lwd=2)
    dev.off()

    # Adding gene symbol and placing it in the front for no and normal shrinkage matricies
    ids.type  <- gene_ID_database_name
    ids <- rownames(res)
    idsN <- rownames(resNorm)
    
    res['gene_id'] <- rownames(res)
    res$gene_name <- as.vector(get.symbolIDsDm(ids,ids.type))
    res_sym_front <- as.data.frame(res) %>% dplyr::select(gene_name, gene_id, everything())    
    
    resNorm['gene_id'] <- rownames(resNorm)
    resNorm$gene_name <- as.vector(get.symbolIDsDm(idsN,ids.type))
    resN_sym_front <- as.data.frame(resNorm) %>% dplyr::select(gene_name, gene_id, everything())    

    # Creating clustermap inputs for no and normal shrinkage
    betasTC <- coef(dds)
    colnames(betasTC)
    topGenes <- which(res_sym_front$padj < PVAL_THRESH, arr.ind = FALSE) # get indicies for all results
    topGenesN <- which(resN_sym_front$padj < PVAL_THRESH, arr.ind = FALSE) # get indicies for all results
    
    # Generating clustermap no shrinkage matrix (NOTE 1,2,3 removed the batch effect columns)
    batch_cols <- seq(1,length(unique(metaData$time))-1)
    mat <- betasTC[topGenes, -(batch_cols)]
    write("batch_cols", stderr())
    write(batch_cols, stderr())
    
    # Generating clustermap normal shrinkage matrix (NOTE 1,2,3 removed the batch effect columns)
    batch_cols <- seq(1,length(unique(metaData$time)))
    matN <- betasTC[topGenesN, -batch_cols]
    
    # Adding gene symbol and placing it in the front for clustermap input matricies (no shrinkage)
    df_mat <- as.data.frame(mat)
    idsM <- rownames(df_mat)
    df_mat['gene_id'] <- idsM
    df_mat$gene_name <- as.vector(get.symbolIDsDm(idsM,ids.type))
    df_mat_sym_front <- as.data.frame(df_mat) %>% dplyr::select(gene_name, gene_id, everything()) 
    
    # Check if NA in gene_name, copy gene_id in its place
    if (NA %in% df_mat_sym_front$gene_name){
        df_mat_sym_front$gene_name <- ifelse(is.na(df_mat_sym_front$gene_name), df_mat_sym_front$gene_id, df_mat_sym_front$gene_name)
    }

    # Adding gene symbol and placing it in the front for clustermap input matricies (normal shrinkage)
    df_matN <- as.data.frame(matN)
    idsMN <- rownames(df_matN)
    df_matN['gene_id'] <- idsMN
    df_matN$gene_name <- as.vector(get.symbolIDsDm(idsMN,ids.type))
    df_matN_sym_front <- as.data.frame(df_matN) %>% dplyr::select(gene_name, gene_id, everything()) 

    # Check if NA in gene_name, copy gene_id in its place
    if (NA %in% df_matN_sym_front$gene_name){
        df_matN_sym_front$gene_name <- ifelse(is.na(df_mat_sym_front$gene_name), df_mat_sym_front$gene_id, df_mat_sym_front$gene_name)
    }

    # Saving sorted (by padj) results
    resSort <- res_sym_front[order(res_sym_front$padj),]
    resSortNormShr <- resN_sym_front[order(resN_sym_front$padj),]
    
    write.csv(as.data.frame(resSort), file=paste(FULL_OUTDIR,paste("deseq2_output_noShrinkage_padj",toString(PVAL_THRESH),'.csv',sep=""), sep="/"), row.names=FALSE, quote=FALSE)
    write.csv(as.data.frame(resSortNormShr), file=paste(FULL_OUTDIR,paste("deseq2_output_normalShrinkage_padj",toString(PVAL_THRESH),'.csv',sep=""), sep="/"), row.names=FALSE, quote=FALSE)
    write.csv(na.omit(df_mat_sym_front),paste(FULL_OUTDIR,paste('deseq2_noShrinkage_clustermapInput_padj',toString(PVAL_THRESH),'.csv',sep=""),sep='/'), row.names=FALSE, quote=FALSE) # clustermap input
    write.csv(na.omit(df_matN_sym_front),paste(FULL_OUTDIR,paste('deseq2_normalShrinkage_clustermapInput_padj',toString(PVAL_THRESH),'.csv',sep=""),sep='/'), row.names=FALSE, quote=FALSE) # clustermap input
}  

if(!interactive()) {
    main()
}