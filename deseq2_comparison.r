# deseq2_comparison.r
# Purpose: run DESeq2 in case vs. control conditions
# Input: 1) metadata file, 2) count matrix file, 3) desired output directory, 4) condition name, and are there 5) batch effects or 6) time components for the experiment set at 7) adj. p-value threshold  
# Output: 1) MA plots (normal and no shrinkage), 2) venn diagram of LRT vs. Wald test (no shrinkage), 3) clustermap differentially expressed gene .csv (normal and no shrinkage), 4) DESeq2 differentially expressed gene output .csvs (normal and no shrinkage) for differentially expressed genes and all genes
# Last Mod: Nov. 19, 2010
# Ashley Conard

#Input arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
      stop("Pass in 9 arguments. Type:
       Rscript DESeq2.r 
            1) /FULL/PATH/TO/METDATA_FILE/ 
            2) /FULL/PATH/TO/COUNT_MATRIX_FILE/ (not normalized and corrected) 
            3) /FULL/PATH/TO/OUTPUTDIR/ (where to create CONDITION_results folder)
            4) CONDITION (e.g. insulin_stim) 
            5) BATCH_EFFECT (0 or 1) 
            6) TIME_COURSE (0 or 1) 
            7) ADJP_THRESH
            8) ORGANISM (dme, hsa, mmu)
            9) CONTROL NAME (from metadata file)")
} else if (length(args) == 9) {
    cat("Passed in:", args,"\n")
} else{
    stop("Pass in 9 arguments. Type:
       Rscript DESeq2.r 
            1) /FULL/PATH/TO/METDATA_FILE/ 
            2) /FULL/PATH/TO/COUNT_MATRIX_FILE/ (not normalized and corrected) 
            3) /FULL/PATH/TO/OUTPUTDIR/ (where to create CONDITION_results folder)
            4) CONDITION (e.g. insulin_stim) 
            5) BATCH_EFFECT (0 or 1) 
            6) TIME_COURSE (0 or 1)             
            7) ADJP_THRESH
            8) ORGANISM (dme, hsa, mmu)
            9) CONTROL NAME (from metadata file)")
}

METDATA <- args[1] 
COUNTDATA <- args[2]
OUTPUTDIR <- args[3] 
CONDITION <- args[4]
BATCH_EFFECT <- as.numeric(args[5]) # 0 or 1, where 1 sets the reduced model to calculate the effect of time, controlling for the effect of batch (in reduced model)
TC <- as.numeric(args[6]) # 0 or 1, if time course then 1
PVAL_THRESH <- as.numeric(args[7])
ORGANISM <- args[8]
CONTROL <- args[9]

# Libraries
library("DESeq2")# make sure 1.22
library(dplyr)
library(annotate)
library(VennDiagram) # to compare Wald and LRT results
source("./clampvsmsl2/geneID_converter.r", local=TRUE)
library(data.table)

# ASSIGNING ORGANISM LIBRARY
if(ORGANISM=="dme"){
    ORG_DB="org.Dm.eg.db"
} else if(ORGANISM=="hsa"){
    ORG_DB="org.Hs.eg.db"
}else if(ORGANISM=="mmu"){
    ORG_DB="org.Mm.eg.db" 
} else{
    stop("Please enter dme (Drosophila melanogaster), hsa (Homo sapiens), or mmu (Mus musculus)")
}
library(ORG_DB, character.only = TRUE) # organism database library
    
# Main function
main <- function(){
    
    # Set database to work with
    gene_ID_database <- as.data.frame(ORG_DB)
    if(ORGANISM=="dme"){
        gene_ID_database_name <- "flybase"
    } else if(ORGANISM!="dme"){
        gene_ID_database_name <- "NA"
    } else{
        stop("Please enter dme (Drosophila melanogaster), hsa (Homo sapiens), or mmu (Mus musculus)")
    }
        
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
    
    
    # Create DESeq2 object and run DESeq2
    # Examine influence of condition (normally what DESeq2 is used for - case vs. control, no time)
    if(BATCH_EFFECT & !TC & CASEvsCONT){ # batch effect and condition
        write("Batch effect and case vs. control", stderr())
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ batch + condition) 
        
        # Set reference level for case control comparison
        dds$condition <- relevel(dds$condition, ref = CONTROL)

        ddsW <- DESeq(dds, test="Wald", betaPrior=TRUE)
        ddsLRT <- DESeq(dds, test="LRT", reduced= ~ batch)
    
    # Examine influence of condition without batch effects (normally what DESeq2 is used for - case vs. control, no time)
    } else if(!BATCH_EFFECT & !TC & CASEvsCONT){ 
        write("Case vs. control", stderr())
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ condition)
        
        # Set reference level for case control comparison
        dds$condition <- relevel(dds$condition, ref = CONTROL)
        
        ddsW <- DESeq(dds, test="Wald", betaPrior=TRUE)
        ddsLRT <- DESeq(dds, test="LRT", reduced= ~1)
        
    # Examine influence of batch effect
    } else if(BATCH_EFFECT & !TC & !CASEvsCONT){
        write("Batch effect", stderr())
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ batch)
        
        # Set reference level for case control comparison
        dds$condition <- relevel(dds$condition, ref = CONTROL)
        
        ddsW <- DESeq(dds, test="Wald", betaPrior=TRUE)
        ddsLRT <- DESeq(dds, test="LRT", reduced= ~1)
    
    # Error
    } else{
        write("ERROR: TIMEOR processes datasets with columns 'time', 'batch', and 'condition'.", stderr())
    }

    # Filtering based on adjusted p-value (i.e. padj)
    res_no_padjW <- results(ddsW)
    resW <- res_no_padjW[which(res_no_padjW$padj < PVAL_THRESH),]
    
    res_no_padjLRT <- results(ddsLRT)
    resLRT <- res_no_padjLRT[which(res_no_padjLRT$padj < PVAL_THRESH),]
    
    # Performing normal shrinkage transformation
    # lfcShrink() should be used downstream of DESeq() with betaPrior=FALSE (the default)
    #resNorm_no_padjW <- lfcShrink(ddsW, coef=3, type="normal") # coef=3 is treatment_odor_vs_etoh  
    #resNormW <- resNorm_no_padjW[which(resNorm_no_padjW$padj < PVAL_THRESH),]
    
    resNorm_no_padjLRT <- lfcShrink(ddsLRT, coef=2, type="normal") # coef=3 is treatment_odor_vs_etoh  
    resNormLRT <- resNorm_no_padjLRT[which(resNorm_no_padjLRT$padj < PVAL_THRESH),]

    # Comparing LRT and Wald test with Venn Diagram
    # Reference: https://cran.r-project.org/web/packages/nVennR/vignettes/nVennR.html
    #            https://thenode.biologists.com/venn-euler-upset-visualize-overlaps-in-datasets/education/
    dfresW <- cbind(IDName = rownames(resW), resW)
    rownames(dfresW) <- 1:nrow(dfresW)
    
    dfresLRT <- cbind(IDName = rownames(resLRT), resLRT)
    rownames(dfresLRT) <- 1:nrow(dfresLRT)
    
    venn.diagram( x = list(as.list(dfresW[["IDName"]]),  as.list(dfresLRT[["IDName"]])), category.names = c("W_DEG" , "LRT_DEG"),filename = paste(FULL_OUTDIR,'venn_wald_LRT.png', sep="/"), output=TRUE) 
    
    # Saving MA plots before and after shrinkage
    pdf(paste(FULL_OUTDIR,paste('ma_plot_noShrinkage_padj_Wald',toString(PVAL_THRESH),'.pdf',sep=""),sep="/"))
    plotMA(resW, ylim=c(-4,4), cex=.8)
    abline(h=c(-1,1), col="dodgerblue", lwd=2)
    dev.off()
#     pdf(paste(FULL_OUTDIR,paste('ma_plot_normalShrinkage_padj_Wald',toString(PVAL_THRESH),'.pdf',sep=""), sep="/"))
#     plotMA(resNormW, ylim=c(-4,4), cex=.8)
#     abline(h=c(-1,1), col="dodgerblue", lwd=2)
#     dev.off()
    pdf(paste(FULL_OUTDIR,paste('ma_plot_noShrinkage_padj_LRT',toString(PVAL_THRESH),'.pdf',sep=""),sep="/"))
    plotMA(resLRT, ylim=c(-4,4), cex=.8)
    abline(h=c(-1,1), col="dodgerblue", lwd=2)
    dev.off()
    pdf(paste(FULL_OUTDIR,paste('ma_plot_normalShrinkage_padj_LRT',toString(PVAL_THRESH),'.pdf',sep=""), sep="/"))
    plotMA(resNormLRT, ylim=c(-4,4), cex=.8)
    abline(h=c(-1,1), col="dodgerblue", lwd=2)
    dev.off()

    # Adding gene symbol and placing it in the front for no and normal shrinkage matricies
    ids.type  <- gene_ID_database_name
    idsW <- rownames(resW)
    idsLRT <- rownames(resLRT)
    #idsN <- rownames(resNormW)
    idsNLRT <- rownames(resNormLRT)
    
    resW['gene_id'] <- rownames(resW)
    resW$gene_name <- as.vector(get.symbolIDsDm(idsW,ids.type))
    res_sym_frontW <- as.data.frame(resW) %>% dplyr::select(gene_name, gene_id, everything())    
    
    #resNormW['gene_id'] <- rownames(resNormW)
    #resNormW$gene_name <- as.vector(get.symbolIDsDm(idsNW,ids.type))
    #resN_sym_frontW <- as.data.frame(resNormW) %>% dplyr::select(gene_name, gene_id, everything()) 
    
    resLRT['gene_id'] <- rownames(resLRT)
    resLRT$gene_name <- as.vector(get.symbolIDsDm(idsLRT,ids.type))
    res_sym_frontLRT <- as.data.frame(resLRT) %>% dplyr::select(gene_name, gene_id, everything())    
    
    resNormLRT['gene_id'] <- rownames(resNormLRT)
    resNormLRT$gene_name <- as.vector(get.symbolIDsDm(idsNLRT,ids.type))
    resN_sym_frontLRT <- as.data.frame(resNormLRT) %>% dplyr::select(gene_name, gene_id, everything()) 

    # Creating clustermap inputs for no and normal shrinkage
    betasTCW <- coef(ddsW)
    colnames(betasTCW)
    topGenesW <- which(res_sym_frontW$padj < PVAL_THRESH, arr.ind = FALSE) # get indicies for all results
    #topGenesNW <- which(resN_sym_frontW$padj < PVAL_THRESH, arr.ind = FALSE) # get indicies for all results
    
    betasTCLRT <- coef(ddsLRT)
    colnames(betasTCLRT)
    topGenesLRT <- which(res_sym_frontLRT$padj < PVAL_THRESH, arr.ind = FALSE) # get indicies for all results
    topGenesNLRT <- which(resN_sym_frontLRT$padj < PVAL_THRESH, arr.ind = FALSE) # get indicies for all results
    
    # Generating clustermap no shrinkage matrix (NOTE 1,2,3 removed the batch effect columns)
    batch_colsW <- seq(1,length(unique(metaData$time))-1)
    matW <- betasTCW[topGenesW, -(batch_colsW)]
    write("batch_colsW", stderr())
    write(batch_colsW, stderr())
    batch_colsLRT <- seq(1,length(unique(metaData$time))-1)
    matLRT <- betasTCLRT[topGenesLRT, -(batch_colsLRT)]
    write("batch_colsLRT", stderr())
    write(batch_colsLRT, stderr())
    
    # Generating clustermap normal shrinkage matrix (NOTE 1,2,3 removed the batch effect columns)
    #batch_colsNW <- seq(1,length(unique(metaData$time)))
    #matNW <- betasTCW[topGenesNW, -batch_cols]
    batch_colsNLRT <- seq(1,length(unique(metaData$time)))
    matNLRT <- betasTCLRT[topGenesNLRT, -batch_colsNLRT]
    
    # Adding gene symbol and placing it in the front for clustermap input matricies (no shrinkage)
    df_matW <- as.data.frame(matW)
    idsMW <- rownames(df_matW)
    df_matW['gene_id'] <- idsMW
    df_matW$gene_name <- as.vector(get.symbolIDsDm(idsMW,ids.type))
    df_mat_sym_frontW <- as.data.frame(df_matW) %>% dplyr::select(gene_name, gene_id, everything()) 
    
    df_matLRT <- as.data.frame(matLRT)
    idsMLRT <- rownames(df_matLRT)
    df_matLRT['gene_id'] <- idsMLRT
    df_matLRT$gene_name <- as.vector(get.symbolIDsDm(idsMLRT,ids.type))
    df_mat_sym_frontLRT <- as.data.frame(df_matLRT) %>% dplyr::select(gene_name, gene_id, everything()) 
    
    # Check if NA in gene_name, copy gene_id in its place
    if (NA %in% df_mat_sym_frontW$gene_name){
        df_mat_sym_frontW$gene_name <- ifelse(is.na(df_mat_sym_frontW$gene_name), df_mat_sym_frontW$gene_id, df_mat_sym_frontW$gene_name)
    }
    
    if (NA %in% df_mat_sym_frontLRT$gene_name){
        df_mat_sym_frontLRT$gene_name <- ifelse(is.na(df_mat_sym_frontLRT$gene_name), df_mat_sym_frontLRT$gene_id, df_mat_sym_frontLRT$gene_name)
    }

    # Adding gene symbol and placing it in the front for clustermap input matricies (normal shrinkage)
#     df_matNW <- as.data.frame(matNW)
#     idsMNW <- rownames(df_matNW)
#     df_matNW['gene_id'] <- idsMNW
#     df_matNW$gene_name <- as.vector(get.symbolIDsDm(idsMNW,ids.type))
#     df_matN_sym_frontW <- as.data.frame(df_matNW) %>% dplyr::select(gene_name, gene_id, everything()) 
    
    df_matNLRT <- as.data.frame(matNLRT)
    idsMNLRT <- rownames(df_matNLRT)
    df_matNLRT['gene_id'] <- idsMNLRT
    df_matNLRT$gene_name <- as.vector(get.symbolIDsDm(idsMNLRT,ids.type))
    df_matN_sym_frontLRT <- as.data.frame(df_matNLRT) %>% dplyr::select(gene_name, gene_id, everything()) 

    # Check if NA in gene_name, copy gene_id in its place
#     if (NA %in% df_matN_sym_frontW$gene_name){
#         df_matN_sym_frontW$gene_name <- ifelse(is.na(df_mat_sym_frontW$gene_name), df_mat_sym_frontW$gene_id, df_mat_sym_frontW$gene_name)
#     }
    if (NA %in% df_matN_sym_frontLRT$gene_name){
        df_matN_sym_frontLRT$gene_name <- ifelse(is.na(df_mat_sym_frontLRT$gene_name), df_mat_sym_frontLRT$gene_id, df_mat_sym_frontLRT$gene_name)
    }

    # Saving sorted (by padj) results
    resSortW <- res_sym_frontW[order(res_sym_frontW$padj),]
    res_no_padjW_sort <- res_no_padjW[order(res_no_padjW$padj),]
    #resSortNormShrW <- resN_sym_frontW[order(resN_sym_frontW$padj),]
    
    resSortLRT <- res_sym_frontLRT[order(res_sym_frontLRT$padj),]
    res_no_padjLRT_sort <- res_no_padjLRT[order(res_no_padjLRT$padj),]
    resSortNormShrLRT <- resN_sym_frontLRT[order(resN_sym_frontLRT$padj),]
    
    write.csv(as.data.frame(resSortW), file=paste(FULL_OUTDIR,paste("deseq2_output_noShrinkage_padj_W",toString(PVAL_THRESH),'.csv',sep=""), sep="/"), row.names=FALSE, quote=FALSE)
    write.csv(as.data.frame(res_no_padjW_sort), file=paste(FULL_OUTDIR,paste("deseq2_output_noShrinkage_padj_W_allGenes",toString(PVAL_THRESH),'.csv',sep=""), sep="/"), row.names=FALSE, quote=FALSE)
    #write.csv(as.data.frame(resSortNormShrW), file=paste(FULL_OUTDIR,paste("deseq2_output_normalShrinkage_padj_W",toString(PVAL_THRESH),'.csv',sep=""), sep="/"), row.names=FALSE, quote=FALSE)
    write.csv(na.omit(df_mat_sym_frontW),paste(FULL_OUTDIR,paste('deseq2_noShrinkage_clustermapInput_padj_W',toString(PVAL_THRESH),'.csv',sep=""),sep='/'), row.names=FALSE, quote=FALSE) # clustermap input
    #write.csv(na.omit(df_matN_sym_frontW),paste(FULL_OUTDIR,paste('deseq2_normalShrinkage_clustermapInput_padj_W',toString(PVAL_THRESH),'.csv',sep=""),sep='/'), row.names=FALSE, quote=FALSE) # clustermap input

    write.csv(as.data.frame(resSortLRT), file=paste(FULL_OUTDIR,paste("deseq2_output_noShrinkage_padj_LRT",toString(PVAL_THRESH),'.csv',sep=""), sep="/"), row.names=FALSE, quote=FALSE)
    write.csv(as.data.frame(res_no_padjLRT_sort), file=paste(FULL_OUTDIR,paste("deseq2_output_noShrinkage_padj_LRT_allGenes",toString(PVAL_THRESH),'.csv',sep=""), sep="/"), row.names=FALSE, quote=FALSE)
    write.csv(as.data.frame(resSortNormShrLRT), file=paste(FULL_OUTDIR,paste("deseq2_output_normalShrinkage_padj_LRT",toString(PVAL_THRESH),'.csv',sep=""), sep="/"), row.names=FALSE, quote=FALSE)
    write.csv(na.omit(df_mat_sym_frontLRT),paste(FULL_OUTDIR,paste('deseq2_noShrinkage_clustermapInput_padj_LRT',toString(PVAL_THRESH),'.csv',sep=""),sep='/'), row.names=FALSE, quote=FALSE) # clustermap input
    write.csv(na.omit(df_matN_sym_frontLRT),paste(FULL_OUTDIR,paste('deseq2_normalShrinkage_clustermapInput_padj_LRT',toString(PVAL_THRESH),'.csv',sep=""),sep='/'), row.names=FALSE, quote=FALSE) # clustermap input
} 

if(!interactive()) {
    main()
}