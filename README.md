# XvsY
## RNA-seq Analysis Pipeline for data across sex and time
CLAMP-null vs CLAMP RNAi vs MSL2 RNAi Effects on Brain Development

![illustration of pipline](https://github.com/ashleymaeconard/XvsY/blob/develop2/detailed_workflow.PNG)

## The Snakemake Pipeline

Snakemake is a workflow management system that was used to create a reproducible and scalable pipeline for RNA-seq analyses, named XvsY, that the user can interact with through the command line. This includes tools to determine the differentially expressed genes based on read count data using DESeq2, find intersections between  gene sets using Intervene, generate global X vs autosome violin plots and comparative fold change box plots, conduct gene ontology analysis using clusterProfiler, and perform motif analysis using MEME and FIMO from the MEMESuite.

To create the environment required to run XvsY, use the following command from the terminal:
```bash
conda env create -f xvsy_environment.yml
```

Also due to large file size, you will need to unzip the dm6.fa and genes.gtf files from the genomes directory with this command:
```bash
unzip genomes/dm6.zip -d genomes/
rm -rf genomes/dm6.zip
unzip genomes/genesgtf.zip -d genomes/
rm -rf genomes/genesgtf.zip
```
also check the gtf file for the genes from way back when

Now to run XvsY with Snakemake use the following commands from the terminal:
```bash
snakemake -R --until RULE --cores 1 --config outdir='/PATH/TO/OUTPUT/DIRECTORY/' ...
```

add rest of config flags
dag visualizations

After the --until flag, you can select any rule from the Snakefile, and the pipeline will run all the rules (with all their scripts) needed to create the inputs for the rule that you would like to run until.

The following is the list of rules in the Snakefile for XvsY:
* run_deseq2 - this rule runs DESeq2 twice to create the sets of differentially expressed genes for two different experiments (eg. RNAi #1 vs RNAi #2 *or* Timepoint #1 vs Timepoint #2 *or* Drug #1 vs Drug #2)
  * config flags:
    * metadata='readcounts/METADATA.CSV' (ex. 'readcounts/metadata_allSamples_full.csv')
    * counts='readcounts/READCOUNTS.CSV' (ex. 'readcounts/FB_gene_ByCondition_countTable_full.csv')
    * condition='CONDITION1' (ex. 'cRNAi_e')
    * control='CONTROL1' (ex. 'elavGRFP_e')
    * condition2='CONDITION2' (ex. 'mRNAi_e')
    * control2='CONTROL2' (ex. 'elavGRFP_e')
    * batch_effect=0 or 1
    * time_course=0 or 1
    * pval_threshold=0.05
    * organism='dme'
    * stat_test='wald' or 'lrt'
    * read_threshold=3
* get_ids - this rule uses the outputs of DESeq2 to generate text files with the gene IDs for the sets of differentially expressed genes
* find_intersections - this rule uses the gene ID text files to find the intersections for the sets of differentially expressed genes
* get_gene - this rule converts the gene ID .txt files to .bed files (gets chromosome, start and end for each gene)
* make_boxplots
* make_XAboxplots
* global_boxplots
* go_analysis
* go_summary
* meme_suite_prep
* run_meme
* run_fimo

<p align="center">
    <img src="https://github.com/ashleymaeconard/XvsY/blob/develop2/snakemake_dag.PNG" width="500" alt="Snakemake Pipeline DAG">
</p>
<!-- ![Snakemake Pipeline DAG](https://github.com/ashleymaeconard/XvsY/blob/develop2/snakemake_dag.PNG) -->

ask ashley if there is a way to create a tree connecting these rules so users can see what rule is required for what
also add in scripts to go from fastq raw data files to read counts?

After the --config flag, set outdir equal to the directory containing the Snakefile, the scripts, and your data, as this is the directory in which XvsY will
generate all of its outputs (*outdir* is the base directory).


## Intersections (BLUE)
First to determine the differentially expressed genes from the RNA-seq read count data, DESeq2 is used in the rule **run_deseq2**.

Using differential expression results from DESeq2, gene sets can be overlapped using Intervene by first running the rule **get_ids** and then **find_intersections**. The rule **get_ids** will generate a text file with the gene IDs (ex. the FlyBase ID for Drosophila genes) from the CSV file outputted by DESeq2. The rule **find_intersections** will then run Intervene on these text files to generate a bar plot to visualize the size of the overlaps as well as text files containing the gene IDs in each overlap.

Next, the rule **get_gene** will generate BED files for the genes in each of these intersections by getting the chromosome number, start site, end site, and gene name.

These scripts will generate Venn Diagrams and Bar Plots to visualize the intersections and they will also generate BED files for the gene sets in each overlap.

Additionally, these intersection scripts rely on individual python scripts:
```bash
python csv_to_id.py -f CSV
```
(to generate text file list of FlyBase IDs from differential expression results CSV file, used in *all_csv_to_id.sh*)
```bash
python gtf_to_bed2.py -f /PATH/TO/genes.gtf -s /PATH/TO/SAVE/genes_gtf.bed
```
(to create BED file from GTF to get gene names, chromosome, start and end, used in *all_id_to_gene.sh*)
```bash
python id_to_gene3.py -t /PATH/TO/TXT_FILE -b /PATH/TO/genes_gtf.bed
```
(to create BED files to Get genes from FlyBase IDs in unique sets from Intervene, used in *all_id_to_gene.sh*)

## Fold Change Comparison between gene sets (GREEN)
To compare the fold changes between two or more of these gene sets in BED file format (from the intersections), use the rule **make_boxplots**.

### X vs A Plots
To generate global X vs A boxplots for a given gene set in BED file format (from the intersections), use the rule **make_XAboxplots**.

To generate violin plots comparing the X vs A fold change between two gene groups in BED file format, use the rule **global_boxplots**.

These rules call the scripts, allBoxPlots.py, xaBoxPlotsfinal.py, and global_XAv5.py, which use the Seaborn package in Python to generate plots to visualize these comparisons.

## GO Analysis (YELLOW)
To run gene ontology (GO) analysis to determine the enriched biological processes (BP), molecular function (MF), and cellular component (CC) of given gene set as a BED file for all BED files in directory, use the rule **go_analysis**.

This rule calls the shell script, ./runGOall.sh, which uses: bed_to_geneListcsv.py and clusterProfiler.r.

To create a summary tables (one for BP, one for MF, and on for CC) containing the top GO term for each of these categories from each gene cluster, use the rule **go_summary**.

## MEME and FIMO Analysis (RED)
For running MEME and/or FIMO for motif analysis, run the rule **meme_suite_prep** to generate the fastq files for each of the gene clusters that is required to run both MEME and FIMO. This rule uses the shell script memesuitePrep.sh and the Python script meme_suite_prep_indiv_cluster.py.

To run MEME analysis to find the top three de novo motifs in each of the gene sets, use the rule **run_meme**. This calls the script runMeme.sh.

To run FIMO analysis to find the presence of CLAMP motifs in each of the gene sets and generate a summary of the FIMO results as a table showing the percentage of a each motif, use the rule **run_fimo**. This calls the script runFimo.sh and fimo_summary.py. NOTE: to run FIMO, a position weight matrix (PWM) is required for your transcription factor or protein of interest.