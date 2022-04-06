<p align="center" width="100%">
  <img src="https://github.com/ashleymaeconard/XvsY/blob/dev_manu_prep/xvsy_logo.PNG" height="100" alt="XvsY Logo">
</p>

## Identify shared and distinct molecular signatures across multiple experiments and time
### Focus on multiple expression and protein-DNA interaction data sets over time
<!---CLAMP-null vs CLAMP RNAi vs MSL2 RNAi Effects on Brain Development in Drosophila-->

![illustration of pipline](https://github.com/ashleymaeconard/XvsY/blob/dev_manu_prep/detailed_workflow.PNG)

**Figure 1**: XvsY enables the user to compare multi-dimensional experiments efficiently and report downstream analyses. First (pink box), given all or selected contexts, the user can choose to run through differential expression analysis and overlap each experiment’s differentially expressed genes (DEGs) to identify distinct and shared features. Second (yellow box), the user can choose to characterize differences, followed by determining mechanism through gene ontology (green box). Lastly (purple box), the user can choose to investigate the upstream gene group targets through motif analysis. The right side illustrates several tool outputs from our proof-of-concept study. Importantly, this analysis can be performed on any experiments that have multiple conditions, and the process can be initiated from anywhere in the pipeline. For example, should the user have differential protein abundances, the user could begin with the intersection step in the pink box. Black boxes denote processing steps. Dotted lines denote an analysis step with output for the user. TF: transcription factor.

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

The `readcounts/` directory contains our proof-of-concept data information. Specifically, 1) the metadata file with four columns: experiment ID, batch, condition name, time point, and 2) the raw read-count data obtained after running HISAT2 and HTSeq. Each column corresponds to a row in the metadata file and each row corresponds to a gene ID.

## Run steps
Run XvsY fully through using the following commands in terminal:
```bash
$snakemake -R --until <RULE> --cores 1 --config outdir='/PATH/TO/OUTPUT/DIRECTORY/' 
```

After the `--until` flag, you can select any rule from the Snakefile, and the pipeline will run through until that specified rule. See [rules](#rules) for details.

## Rules
Here are the list of XvsY rules in the Snakefile workflow, and their corresponding `config` flags. Note acronyms and definitions. *DEGs*: differentially expressed genes; *GO*: gene ontology; *context*: experiments, replicates and conditions in those experiments, and any time points or other conditions such as sex in those conditions; *group*: the resulting DEG sets after finding intersections (i.e. overlaps); *distinct* DEGs: the set of DEGs beloning only to that context; *shared* DEGs: the set of DEGs beloning to two or more contexts. See Figure 1 above for pictoral  representation.

* `run_deseq2` 
  * runs DESeq2 for each set of experiments independently, so to compare input contexts (e.g. experiment 1's case vs. control, then experiment 2's case vs. control)
  * `config` flags:
    * `num_comps`: 2 or 3 
        * number of experiments/time points to be compared
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
  * uses each .bed file to compare the fold-changes between all pairs of DEG groups and DEG contexts
* `make_XAboxplots`
  * uses each .bed file to compare the fold-changes of the DEGs on the X chromosome vs. the autosomes between contexts.
* `global_boxplots` 
  * uses DESeq2 outputs to generate global context chromosome notch plots depicting fold-changes across all chromosomes
* `go_analysis`
  * uses each .bed file to run gene ontology (GO) analysis for each group
  * `config` flags:
    * `sep_tps`: 0 or 1 
      * 1 to run GO for each time point separately, 0 otherwise
      * stands for 'separate time points'.
    * `organism`: 'dme'
    * `go_pval`: 0.05
      * any continuous value between (0-1]
* `go_summary` 
  * collates information across all groups GO analyses to create a summary table per category (biological process, molecular function, cellular component).
* `meme_suite_prep` 
  * uses .bed files to generate .fastq files for input into MEME and/or FIMO
  * `config` flags:
    * `reform_genes`: <PATH/TO/XvsY/genomes/REFORMATTED_GENES.CSV>
        * example: 'genomes/reformatted_genes_gtf.csv'
    * `chrom_fa`: <PATH/TO/genomes/ORG.FA> 
        * example: 'genomes/dm6.fa'
    * `tss_only`: 1 or 0 
        * 1 to run +-1kb from transcription start site, 0 to consider gene body region.
    *  organism='dme'
* `run_meme`
    * uses .fastq files to run MEME and return the top 3 *de novo* motifs
* `run_fimo` 
    * uses 1) .fastq files and 2) input position weight matrix (PWM) to run FIMO and find TF binding motifs across groups of DEGs
    * `config` flags:
      * `pwm_path`: <PATH/TO/pwms/PWM_OF_INTEREST.txt> 
        * 'pwms/meme_CLAMP_overlap_GAF.txt'

<p align="center">
    <img src="https://github.com/ashleymaeconard/XvsY/blob/dev_manu_prep/snakemake_dag.PNG" width="500" alt="Snakemake Pipeline DAG">
</p>

**Figure 2**: The XvsY Snakemake workflow can be shown as a directed acyclic graph of XvsY rules to allow the user to start at any rule (rectangle), given that its inputs are present, and run until any connecting rule. The colors are associated with the steps in Figure 1. The paths illustrate the rules to be executed and their associated dependencies.

Once the `--config flag` is set, make `outdir` equal to the directory containing the Snakefile, the scripts, and the data, as this is the directory in which XvsY will generate all of its outputs. That means that `outdir` is the base directory where all results will be stored.

## XvsY example/ demo run
For example, to run the XvsY Snakemake pipeline until the `go_summary` rule for GO analysis, assuming the user is in the cloned XvsY directory, run the following one line command from the terminal:
```bash
$ snakemake -R --until go_summary --cores 1 --config outdir='~/xvsy_test/test/' metadata='readcounts/metadata_allSamples_full.csv' counts='readcounts/FB_gene_ByCondition_countTable_full.csv' condition='cRNAi_e' condition2='mRNAi_e' batch_effect=0 time_course=0 pval_threshold=0.05 organism='dme' control='elavGRFP_e' control2='elavGRFP_e' stat_test='wald' read_threshold=3 gtf_path='genomes/genes.gtf' sep_tps=0 organism='dme' go_pval=0.05
```
Should the user prefer to run individual rules, below are specifications for which rule to run, and when. 

## Identify distinct and shared groups through differential expression and intersections (pink)
Use rule `run_deseq2` to first determine differentially expressed genes (DEGs) from the RNA-seq read count data using DESeq2.

Using those DEG results per context, each context can then be overlapped using Intervene by first running the rule `get_ids` and then `find_intersections`. 
`get_ids` generates a text file with the gene IDs (eg. the FlyBase ID for *Drosophila* genes) from the .csv file that is output by DESeq2. 
`find_intersections` then runs Intervene on these text files to generate an overlap to both visualize the overlaps and create text files of the created groups (i.e. gene IDs in each overlap).

Next, the rule `get_gene` generates .bed files for each group by getting the gene name, chromosome number, start site, and end site, in that order.

## Uncover fold-change differences for distinct and shared groups (yellow)

To generate global context chromosome notch plots for a given context, use the  `make_XAboxplots` rule. NOTE, while the rule says box plots, they are now notch plots. This verbage will be updated in a future version.

Use the `make_boxplots` rule to compare the fold-changes between two or more of these groups. NOTE, while the rule says box plots, they are now violin plots. This verbage  will be updated in a future version.

To generate violin plots comparing the X vs. A fold-change between two gene groups, use the rule `global_boxplots`. NOTE, while the rule says box plots, they are now violin plots. This verbage  will be updated in a future version.

## Characterize mechanism for distinct and shared groups (green)
To run gene ontology (GO) analysis to determine the enriched biological processes (BP), molecular function (MF), and cellular component (CC) of given DEG group or for all groups in a directory, use the rule `go_analysis`.

Then, use the `go_summary` rule to create summary tables for each category (one for BP, MF, and CC separately) containing the top GO term for each group.

## Find motifs and binding regions for distinct and shared groups (purple)
To perform *de novo* motif discovery and motif finding for a given transcription factor (TF) of interest, first use `meme_suite_prep` to generate the .fastq files for each group. These are needed to subsequently run MEME and/or FIMO.

To run MEME analysis which returns the top three *de novo* motifs for each DEG group, use the rule `run_meme`.

To run FIMO analysis which searches across group DEG members to find individual matches to the TF motif provided, use the rule `run_fimo`. This rule requires an input TF position weight matrix (PWM) for the TF of interest. XvsY then generates a summary of the FIMO results as a table showing the percentage of a each motif.

## Choosing between tests
In the second stage, each group is compared intra and inter-context (yellow, Figure 1). Specifically, XvsY outputs: 1) global context chromosome level notch plots with a heatmap significance table highlighting significant expression differences across chromosomes; 2) expression fold-change violin plots for shared group differentially expressed genes; and 3) X vs. autosome violin plots for genes that are upregulated and downregulated within and between group (Figure 2).  For each type of plot, there are three significance measures that are output. For example, in the first global context chromosome level notch plot, the significance of the difference in expression between each pair of chromosomes is computed using 3 tests, the t-test, Kolmogorov-Smirnov (KS) test, and Mann-Whitney (MW) U test (also called Wilcoxon Rank-Sum test). The resulting p-values are displayed in a heatmap with the darker values indicating stronger significance. Due to the issue of multiple comparisons, all p-values are adjusted using the Benjamini-Hochberg (BH) correction method. The user has the choice of deciding which test best fits their hypothesis. 

Specifically, the user is encouraged to examine the most appropriate test statistic or p-value summaries for their analysis when determining the fold-change differences between chromosomes. Briefly, the t-test provides an exact test to determine if the means of two i.i.d. normal distributions (with equal and unknown variances) are the same. When the normality assumption does not hold, a non-parametric alternative such as the MW or KS test will likely have better statistical power. Nevertheless, with differing variances between groups, a t-test could provide better type-1 error control, that is, control of false positives where one would reject a true null hypothesis. The appropriate test should be chosen based on the stated hypothesis.

Non-parametric tests do not typically test for mean differences, so both the MW and KS tests should be used carefully if mean differences are of primary interest to the user. The MW test is a nonparametric test of the null hypothesis that for any randomly selected values *i* and *j* from two populations, the probability *Pr(i>j)=Pr(j>i)*. The Kolmogorov-Smirnov test is a nonparametric test comparing the equality of two one-dimensional probability distributions. The null hypothesis states that the two sets of samples were drawn from the same probability distribution. Although the t-test and MW test are more commonly used in the analysis of biological data, they make assumptions about the distributions of the data and the sample sizes that may not always be justified. For example, while the MW test is sensitive to changes in the median, if the user is interested in  substantial differences in the shape or spread between two distributions, the KS test is more sensitive. For these reasons, XvsY includes several statistical tests and provides guidance regarding which test to use.