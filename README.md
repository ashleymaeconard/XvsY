<p align="center" width="100%">
  <img src="https://github.com/ashleymaeconard/XvsY/blob/dev_manu_prep/xvsy_logo.PNG" height="100" alt="XvsY Logo">
</p>

## Identify shared and distinct molecular signatures across multiple experiments and time
### Focus on multiple expression and protein-DNA interaction data sets over time
<!---CLAMP-null vs CLAMP RNAi vs MSL2 RNAi Effects on Brain Development in Drosophila-->

![illustration of pipline](https://github.com/ashleymaeconard/XvsY/blob/dev_manu_prep/detailed_workflow.PNG)

### **Figure 1**: XvsY enables the user to compare multi-dimensional experiments efficiently and report downstream analyses. First (pink box), given all or selected contexts, the user can choose to run through differential expression analysis and overlap each experiment’s differentially expressed genes (DEGs) to identify distinct and shared features. Second (yellow box), the user can choose to characterize differences, followed by determining mechanism through gene ontology (green box). Lastly (purple box), the user can choose to investigate the upstream gene group targets through motif analysis. The right side illustrates several tool outputs from our proof-of-concept study. Importantly, this analysis can be performed on any experiments that have multiple conditions, and the process can be initiated from anywhere in the pipeline. For example, should the user have differential protein abundances, the user could begin with the intersection step in the pink box. Black boxes denote processing steps. Dotted lines denote an analysis step with output for the user. TF: transcription factor.

## Quick Start: 3 Steps

  1. Activate XvsY environment using the following command:
  ```
  conda activate xvsy_environment
  ```
  2. Locate your data and metadata file directory (<PATH/TO/OUTPUT/DIRECTORY>). The data comprise of a 1) read-count matrix of genes by contexts and 2) associated metadata .csv file. See [Inputs](#inputs-for-xvsy) section for details.

  3. Run this command within `XvsY/` directory to run all the way through XvsY (four boxes in Figure 1 above):
  ```
  snakemake --R --until <RULE> --cores 1 --config outdir=<PATH/TO/OUTPUT/DIRECTORY>
  ```

## Paper

Conard, A. M., Nathoo, I., Tsiarli, M. A., Larschan, E. N. (2022). XvsY: a tool to identify shared and distinct molecular signatures across multiple experiments and time.


## Overview

We present XvsY, the first method to comprehensively analyze complex multi-dimensional experiments to identify distinct and shared signatures and associated mechanisms. This is a Snakemake pipeline command-line interface (CLI) tool that includes four stages: 1) identifying distinct and shared features through differential expression and/or intersections to form significant gene groups; 2) characterizing fold-change differences between experimental groups; 3) uncovering gene function; 4) discovering DNA sequence motifs for all groups (Figure 1). Should the user choose to run through the entire analysis, follow the [Quick Start](#quick-start-3-steps).


## Installation

Assuming you have [Anaconda](https://docs.anaconda.com/anaconda/install/index.html) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed, please follow the commands below:

1. Clone the Github repository either locally or on a remote server.
2. Navigate to the XvsY directory in terminal:
  ```bash
  $ cd /CLONED/REPO/PATH/XvsY/
  ```
3. Activate your base conda environment (example using Miniconda):
  ```bash
  $ source <DIR>/bin/activate
  ```
4. Type the following command in the current directory
  ```bash
  $ conda env create -f xvsy_environment.yml
  ```
5. Activate your XvsY environment:
  ```bash
 $ conda activate xvsy_environment
  ```
6. Unzip the `dm6.fa` and `genes.gtf` files within the genomes directory, from the `XvsY/` directory.
```bash
$ unzip genomes/dm6.zip -d genomes/
$ rm -rf genomes/dm6.zip
$ unzip genomes/genesgtf.zip -d genomes/
$ rm -rf genomes/genesgtf.zip
```


## Inputs for XvsY

XvsY requires a read-count matrix with the gene IDs as rows and context IDs (i.e. experiments, replicates and conditions in those experiments, and any time points or other conditions such as sex in those conditions) IDs as columns. Further, a metadata file is needed to specify the rows as the read-count matrix's context IDs, and the colums as context ID, condition, batch, and time point. Each row of the metadata file corresponds to a column in the read-count matrix.

## Run XvsY

### Framework

XvsY is an accessible command-line interface (CLI) enabling the user to accessibly and efficiently compare multi-dimensional experiments, with a focus on RNA-seq. This tool is suitable for the analysis of any experimental contexts that produce differentiable data, such as RNA interference, CRISPR, or therapeutic drug testing.
We adopted Snakemake as a workflow management system to create a reproducible and scalable pipeline. 

### Directories

The `scripts/` directory contains all the Python, R, and bash scripts required to run XvsY. The genomes directory contains the Drosophila genes GTF file, the chromosome FA file, and the reformatted genes file (from the GTF). The PWMs directory contains the PWMs required for FIMO analysis (if desired).

The `readcounts/` directory contains our proof-of-concept data information. Specifically, 1) the metadata file with four columns: experiment ID, batch, condition name, timepoint, and 2) the raw read-count data obtained after running HISAT2 and HTSeq. Each column corresponds to a row in the metadata file and each row corresponds to a gene ID.

## Run steps
Run XvsY fully through using the following commands in terminal:
```bash
$snakemake -R --until <RULE> --cores 1 --config outdir='/PATH/TO/OUTPUT/DIRECTORY/' 
```

After the `--until` flag, you can select any rule from the Snakefile, and the pipeline will run through until that specified rule. See [rules](#rules) for details.

## Rules
Here are the list of XvsY rules in the Snakefile workflow, and their corresponding `config` flags. Note acronyms and definitions. *DEGs*: differentially expressed genes; *context*: experiments, replicates and conditions in those experiments, and any time points or other conditions such as sex in those conditions; *group*: the resulting DEG sets after finding intersections (i.e. overlaps); *distinct* DEGs: the set of DEGs beloning only to that context; *shared* DEGs: the set of DEGs beloning to two or more contexts. See Figure 1 above for pictoral  representation.

* `run_deseq2` 
  * runs DESeq2 for each set of experiments independently, so to compare input contexts (e.g. experiment 1's case vs. control, then experiment 2's case vs. control)
  * `config` flags:
    * `num_comps`: 2 or 3 
        * number of experiments/timepoints to be compared
    * `metadata`: <PATH/TO/METADATA>
        * example 'readcounts/metadata_allSamples_full.csv'
    *  `counts`: <PATH/TO/READ-COUNTS>
        * example: 'readcounts/FB_gene_ByCondition_countTable_full.csv')
    * `condition`: 'CONDITION1'
        * example: 'cRNAi_e'
    * `control`: 'CONTROL1'
        * example: 'elavGRFP_e'
    * `condition2`: 'CONDITION2' 
        * example: 'mRNAi_e'
    * `control2`: 'CONTROL2' 
        * example: 'elavGRFP_e'
    * `batch_effect`: 0 or 1
        * 0 for no, 1 for yes
    * `time_course`: 0 or 1
        * 0 for no, 1 for yes
    * `pval_threshold`: 0.05
        * any continuous value between (0-1]
    * `organism`: 'dme'
    * `stat_test`: 'wald' or 'lrt'
        * wald is Wald Test and lrt is likelihood-ratio test
    * `read_threshold`: 3
        * threshold for minimum number of reads to consider a gene
* `get_ids`
  * uses DESeq2 outputs to generate text files with gene IDs for each context of differentially expressed genes (DEGs)
* `find_intersections` 
  * uses the gene ID text files to find the intersections between each context of differentially expressed genes
* `get_gene`
  * converts gene ID text files into .bed files 
    * .bed file columns:  gene ID, chromosome, start, end 
  * `config` flags:
    * `gtf_path`: <PATH/TO/XvsY/genomes/GENES.GTF> 
      * example: 'genomes/genes.gtf'
      * this is the location of the genome annotation file
* `make_boxplots`
  * uses each .bed file to compare the fold changes between all pairs of DEG groups and DEG contexts
* `make_XAboxplots`
  * uses each .bed file to compare the fold changes of the DEGs located on the X chromosome vs the autosomes between different sets
* `global_boxplots` 
  * uses the outputs of deseq2 to generate global X vs A plots and look at the fold changes across all chromosomes
* go_analysis - this rule uses the .bed files to run gene ontology analysis for each gene set
  * config flags:
    * sep_tps=0 or 1 (set to 1 to run GO for each timepoint separately, 0 otherwise)
    * organism='dme'
    * go_pval=0.05
* go_summary - this rule uses the information from each GO cluster to generate summary tables for each GO category
* meme_suite_prep -  this rule uses the .bed files to generate FASTA files for the genes and other inputs for MEME/FIMO
  * config flags:
    * reform_genes='genomes/REFORMATTED_GENES.CSV' (eg. 'genomes/reformatted_genes_gtf.csv')
    * chrom_fa='genomes/ORG.FA' (eg. 'genomes/dm6.fa')
    * tss_only=1 (set to 1 to run +-1kb from transcription start site, 0 otherwise)
    *  organism='dme'
* run_meme - this rule uses the FASTA files to run MEME and find the top 3 de novo motifs
* run_fimo - this rule uses the FASTA files to run FIMO and find the TF binding motif across all genes based on the PWM
  * config flags:
    * pwm_path='pwms/PWM_OF_INTEREST.txt' (eg. 'pwms/meme_CLAMP_overlap_GAF.txt')

<p align="center">
    <img src="https://github.com/ashleymaeconard/XvsY/blob/dev_manu_prep/snakemake_dag.PNG" width="500" alt="Snakemake Pipeline DAG">
</p>
<!-- ![Snakemake Pipeline DAG](https://github.com/ashleymaeconard/XvsY/blob/dev_manu_prep/snakemake_dag.PNG) -->

<!-- ask ashley if there is a way to create a tree connecting these rules so users can see what rule is required for what -->
<!-- also add in scripts to go from fastq raw data files to read counts? -->

After the --config flag, set outdir equal to the directory containing the Snakefile, the scripts, and your data, as this is the directory in which XvsY will generate all of its outputs (*outdir* is the base directory).

For egample, to run the XvsY Snakemake Pipeline until the go_summary rule, run the following command from the terminal:
```bash
snakemake -R --until go_summary --cores 1 --config outdir='/data/compbio/inathoo/snakemake_test/test10/' metadata='readcounts/metadata_allSamples_full.csv' counts='readcounts/FB_gene_ByCondition_countTable_full.csv' condition='cRNAi_e' condition2='mRNAi_e' batch_effect=0 time_course=0 pval_threshold=0.05 organism='dme' control='elavGRFP_e' control2='elavGRFP_e' stat_test='wald' read_threshold=3 gtf_path='genomes/genes.gtf' sep_tps=0 organism='dme' go_pval=0.05
```

## Intersections (BLUE)
First to determine the differentially expressed genes from the RNA-seq read count data, DESeq2 is used in the rule **run_deseq2**.

Using differential egpression results from DESeq2, gene sets can be overlapped using Intervene by first running the rule **get_ids** and then **find_intersections**. The rule **get_ids** will generate a text file with the gene IDs (eg. the FlyBase ID for Drosophila genes) from the CSV file outputted by DESeq2. The rule **find_intersections** will then run Intervene on these text files to generate a bar plot to visualize the size of the overlaps as well as text files containing the gene IDs in each overlap.

Negt, the rule **get_gene** will generate BED files for the genes in each of these intersections by getting the chromosome number, start site, end site, and gene name.

These scripts will generate Venn Diagrams and Bar Plots to visualize the intersections and they will also generate BED files for the gene sets in each overlap.

Additionally, these intersection scripts rely on individual python scripts:
```bash
python csv_to_id.py -f CSV
```
(to generate text file list of FlyBase IDs from differential egpression results CSV file, used in *all_csv_to_id.sh*)
```bash
python gtf_to_bed.py -f /PATH/TO/genes.gtf -s /PATH/TO/SAVE/genes_gtf.bed
```
(to create BED file from GTF to get gene names, chromosome, start and end, used in *all_id_to_gene.sh*)
```bash
python id_to_gene_map.py -t /PATH/TO/TXT_FILE -b /PATH/TO/genes_gtf.bed
```
(to create BED files to Get genes from FlyBase IDs in unique sets from Intervene, used in *all_id_to_gene.sh*)

## Fold Change Comparison between gene sets (GREEN)
To compare the fold changes between two or more of these gene sets in BED file format (from the intersections), use the rule **make_boxplots**.

### X vs A Plots
To generate global X vs A boxplots for a given gene set in BED file format (from the intersections), use the rule **make_XAboxplots**.

To generate violin plots comparing the X vs A fold change between two gene groups in BED file format, use the rule **global_boxplots**.

These rules call the scripts, allBoxPlots.py, xaBoxPlotsfinal.py, and global_XA.py, which use the Seaborn package in Python to generate plots to visualize these comparisons.

## GO Analysis (YELLOW)
To run gene ontology (GO) analysis to determine the enriched biological processes (BP), molecular function (MF), and cellular component (CC) of given gene set as a BED file for all BED files in directory, use the rule **go_analysis**.

This rule calls the shell script, ./runGOall.sh, which uses: bed_to_geneListcsv.py and clusterProfiler.r.

To create a summary tables (one for BP, one for MF, and on for CC) containing the top GO term for each of these categories from each gene cluster, use the rule **go_summary**.

## MEME and FIMO Analysis (RED)
For running MEME and/or FIMO for motif analysis, run the rule **meme_suite_prep** to generate the fastq files for each of the gene clusters that is required to run both MEME and FIMO. This rule uses the shell script memesuitePrep.sh and the Python script meme_suite_prep_indiv_cluster.py.

To run MEME analysis to find the top three de novo motifs in each of the gene sets, use the rule **run_meme**. This calls the script runMeme.sh.

To run FIMO analysis to find the presence of CLAMP motifs in each of the gene sets and generate a summary of the FIMO results as a table showing the percentage of a each motif, use the rule **run_fimo**. This calls the script runFimo.sh and fimo_summary.py. NOTE: to run FIMO, a position weight matrix (PWM) is required for your transcription factor or protein of interest.
