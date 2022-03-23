# clusterProfiler.r
# Purpose: runs clusterProfiler to find and visualize gene enrichment within each cluster.
# Reference: TIMEOR
# Last mod. 02/22/2022

# CHECKING ARGUMENTS
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Type: /usr/local/bin/Rscript clusterProfiler.r
            1) /PATH/TO/INPUT_OUTPUT_DIR/ (e.g. /timeor/results/primary/<EXPERIM_results/clusters/1/)
            2) EXPERIMENT_NAME (e.g. test_results, write 'test')
            3) SEPARATE_TIMEPOINTS (set to 1 to run GO for each timepoint separately, 0 otherwise)
            4) ORGANISM (dme, hsa, mmu)
            5) ADJ_PVAL (recommend 0.05)", call.=FALSE)
} else if (length(args) == 5) {
    cat("Passed in:", args,"\n")
} else{
    stop("Pass in 1) /PATH/TO/INPUT_OUTPUT_DIR/ (e.g. /timeor/results/primary/<EXPERIM_results/clusters/1/)
                  2) EXPERIMENT_NAME (e.g. test_results, write 'test')
                  3) SEPARATE_TIMEPOINTS (set to 1 to run GO for each timepoint separately, 0 otherwise)
                  4) ORGANISM (dme (Drosophila melanogaster), hsa (Homo sapiens), mmu (Mus musculus)
                  5) ADJ_PVAL (recommend 0.05))")
}

# ASSIGNING INPUT ARGUMENTS
IN_OUTPUT <- args[1] 
EXPERIMENT_NAME <- args[2] # without "_results" added
SEP_TPS <- args[3] # 0 or 1
ORGANISM <- args[4] # dme (Drosophila melanogaster), hsa (Homo sapiens), mmu (Mus musculus)
options(digits=5)
ADJ_PVAL <- as.double(args[5]) # recommend 0.05

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

# LOADING PACKAGES
library(biomaRt)
library(devtools)
#install_github("GuangchuangYu/bitr")
library(bitr)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(ORG_DB, character.only = TRUE) # organism database library
library(cowplot)
library(ggplot2)
library(stringr)
library(scales)

print("Done Loading!")

assess_enrichment <- function(geneENS, type_enr, keyz){

    # GO over-representation (enrichement) test
    print(geneENS)
    ego <- enrichGO(gene        = geneENS,
                    OrgDb         = ORG_DB,
                    keyType       = keyz,
                    ont           = type_enr,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = ADJ_PVAL,
                    qvalueCutoff  = ADJ_PVAL)
    print("im here")
    if(is.null(ego)){
	df <- data.frame()
	return(df)
    }
    ego <- setReadable(ego, OrgDb = ORG_DB)
    return(ego)
}

###define a new dotplot function
## x is the variable of the enrichment results
## width is used to set the width of the labels of y-axis
## top is used to set how many toppest terms will be shown in the figure
dotplot_ylab <- function(x,width,top){
	x.df <- as.data.frame(x)
	if(nrow(x.df)>top){
		x.df <- x.df[1:top,]
		}
	x.df$GeneRatio <- unlist(lapply(as.list(x.df$GeneRatio), function(x) eval(parse(text=x))))
	x.df$Description <- factor(x.df$Description,levels = x.df$Description[order(x.df$Count,decreasing = F)])
	p<-ggplot(x.df, aes(x=GeneRatio, y=Description,color=p.adjust)) + theme_bw() +
	geom_point(aes(size = Count))+scale_color_gradient(low="red", high="blue")+scale_y_discrete(labels=function(x) str_wrap(x,width=width))
}

generate_plots <- function(egoo, gl, ex, type_enrich, outdir, num_gene_cluster){
    # Checking if the enrichment object is empty, and if so, exit function and script
    if(nrow(egoo)>0){

        # Generating subfolder within specific cluster
        subdirec <- file.path(outdir, type_enrich)
        if (!dir.exists(subdirec)){
            dir.create(subdirec)
        } else {
            cat(subdirec, "subdirectory exists.")
        }
        # Generate www folder within each subfolder for dotplot rendering in R Shiny
        wwwsvg_subdirec <- file.path(subdirec, "www")
        if (!dir.exists(wwwsvg_subdirec)){
            dir.create(wwwsvg_subdirec)
        } else {
            cat(wwwsvg_subdirec, "www/ subdirectory for *.svg dotplot exists.")
        }

        # Plotting GOrilla plot
        write.csv(egoo, paste(subdirec, paste(ex,"cluster_", num_gene_cluster, "clustProf",type_enrich, ".csv", sep="_"), sep="/"))
        pdf(paste(paste(subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_dag",paste(type_enrich,".pdf",sep=""), sep="_"), sep="/")), width = 10, height = 10)
        plt_plotGOgraph <- plotGOgraph(egoo)
        print(plt_plotGOgraph)
        dev.off()

        # Plotting dot plot
        png(paste(paste(wwwsvg_subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_dotplot",paste(type_enrich,".png",sep=""), sep="_"), sep="/")))
        plt_dotplot <- dotplot_ylab(egoo,40,10)
        print(plt_dotplot)
        dev.off()

        # Plotting clusterProfiler relationship graph plot
        pdf(paste(paste(subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_emaplot",paste(type_enrich,".pdf",sep=""), sep="_"), sep="/")), width = 10, height = 10)
        plt_emaplot <- tryCatch({
	   emapplot(egoo)
	}, error = function(err) {
	   barplot(egoo, showCategory=20)
	})
        print(plt_emaplot)
        dev.off()

        # Barplot
        pdf(paste(paste(subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_barplot",paste(type_enrich,".pdf",sep=""), sep="_"), sep="/")), width = 10, height = 10)
        plt_bar<- barplot(egoo, showCategory=10)
        print(plt_bar)
        dev.off()

        # Gene concept network
        pdf(paste(paste(subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_conceptplot",paste(type_enrich,".pdf",sep=""), sep="_"), sep="/")), width = 30, height = 15)
        p2 <- cnetplot(egoo, categorySize="pvalue", foldChange=gl)
        plt_cnet <- cnetplot(egoo, foldChange=gl, categorySize="pvalue", colorEdge = TRUE)
        plt_cnet_cir <- cnetplot(egoo, foldChange=gl, circular = TRUE, categorySize="pvalue", colorEdge = TRUE)
        plt_concept <- plot_grid(plt_cnet, plt_cnet_cir, ncol=2)
        print(plt_concept)
        dev.off()

        # Plotting clusterProfiler phylogeny plot (if there is more than 1 enriched GO term)
        if(dim(egoo)[1]>1){
            pdf(paste(paste(subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_goplot",paste(type_enrich,".pdf",sep=""), sep="_"), sep="/")), width = 10, height = 10)
            plt_goplot <- goplot(egoo)
            print(plt_goplot)
            dev.off()
        }
    } else{
        print(paste(type_enrich," is empty", sep=""))
    }
}

main <- function(){
    # Finding only filenames that are *geneList*.csv or *geneName*.csv to get list of genes to calc. enrichment
    if (SEP_TPS == 0){
        dir_file <- Sys.glob(file.path(IN_OUTPUT, "*geneList*.csv"))
    } else {
        dir_file <- Sys.glob(file.path(IN_OUTPUT, "*geneName*.csv"))
    }

    if (! length(dir_file)>1){ # check to make sure only one geneList file per cluster
        geneNumClust=(sub(".*_([^.]+)\\.csv.*", "\\1", dir_file))
        OUTPUT_DIR <-paste(dirname(dir_file),'',sep="/")
        print(dir_file)
        # Generating ranked gene list for gene set enrichment analysis (GSEA)
        f <- read.csv(dir_file, sep=",", header = TRUE)
        geneList<-f[,2]
        names(geneList)<- as.character(f[,1])
        geneList<-sort(geneList, decreasing=TRUE)

        # Reading the gene file for each cluster
        a <- dput(as.character(f[1][,1]))

        # Removing ballgown ID and/or transcript ID if present at end of gene name
        if (grepl("_", a[1])){
            cat("Removing values after and including '_' at end of gene/transcript name\n")
            xs <- gsub("\\-R[A-Z][_|\\>].*","",a)
        } else if (grepl("-R[A-Z]\\>", a[1])) {
            cat("Removing transcript ID at end of gene name\n")
            xs <- gsub("\\-R[A-Z]\\>.*","",a)
        } else{
            xs <- a
        }

        # Converting gene names to ensemble and entrezids
        gene = bitr(geneID=xs, fromType="SYMBOL", toType=c("ENSEMBL","ENTREZID"), OrgDb=ORG_DB)
        if ("FBgn0003996" %in% xs) {
            print("Adding white gene")
            gene <- gene[,c("SYMBOL","ENSEMBL","ENTREZID")]
            bbbbb <- data.frame("w","FBgn0003996","31271")
            names(bbbbb) <- c("SYMBOL","ENSEMBL","ENTREZID")
            gene <- rbind(gene,bbbbb)
        }
        #gene <- gene[,c("SYMBOL","ENSEMBL","ENTREZID")]
        #bbbbb <- data.frame("w","FBgn0003996","31271")
        #names(bbbbb) <- c("SYMBOL","ENSEMBL","ENTREZID")
        #gene <- rbind(gene,bbbbb)
        #gene <- gene[-131,]
        print(dim(gene))
        print(head(gene,150))
        print("DONE WITH BITR")

        # Writing an intermediate file with gene names
        type_enrichment <- c("BP","MF","CC") # Biological Process (BP), Molecular Function (MF), Cellular Component (CC)

        # Plotting GO and pathway enrichment per type (BP, MF, CC)
        pathway_found=0 # flip to 1 when the most enriched pathway is identified
        for(types in type_enrichment){
            print(paste("Processing", geneNumClust, sep=" "))
            print(paste("Determining",types,"enrichment", sep=" "))
	    print(gene)
	    ego <- assess_enrichment(gene$ENSEMBL, types, 'ENSEMBL')
            #ego3 <- assess_GSEA(geneList,type)

            print(paste("Plotting",types,"enrichment for", EXPERIMENT_NAME,sep=" "))
	    if(nrow(ego)==0){
		    (print(paste("Trying enrichment with ENTREZID ID")))
		     ego <- assess_enrichment(gene$ENTREZID, types, 'ENTREZID')
	    }
            generate_plots(ego, geneList, EXPERIMENT_NAME, types, OUTPUT_DIR, geneNumClust)
        }

    }
}

if(!interactive()) {
    main()
}
