<p align="center" width="100%">
  <img src="https://github.com/ashleymaeconard/XvsY/blob/dev_manu_prep/xvsy_logo.PNG" height="100" alt="XvsY Logo">
</p>

## Identify shared and distinct genomic signatures from multiple expression and protein-DNA interaction data sets over time
<!---CLAMP-null vs CLAMP RNAi vs MSL2 RNAi Effects on Brain Development in Drosophila-->

![illustration of pipline](https://github.com/ashleymaeconard/XvsY/blob/dev_manu_prep/detailed_workflow.PNG)

## Quick Start: 3 Steps

## Paper and Citation

## Overview

## Installation

## Run XvsY

## Details

### Snakemake framework

We adopted Snakemake as a workflow management system to create a reproducible and scalable pipeline. that the user can interact with through the command line. This includes tools to determine the differentially egpressed genes based on read count data using DESeq2, find intersections between  gene sets using Intervene, generate global X vs autosome violin plots and comparative fold change box plots, conduct gene ontology analysis using clusterProfiler, and perform motif analysis using MEME and FIMO from the MEMESuite.

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
<!-- also check the gtf file for the genes from way back when -->

The scripts directory contains all the Python, R, and bash scripts required to run XvsY. The genomes directory contains the Drosophila genes GTF file, the chromosome FA file, and the reformatted genes file (from the GTF). The PWMs directory contains the PWMs required for FIMO analysis (if desired).

In regards to the data, the readcounts directory contains the metadata file with four columns: egperiment ID, batch, condition name, timepoint. It also includes the raw read count data obtained after running HISAT2 and HTSeq. This time has a column corresponding to each of the egperiment IDs in the metadata file and a row for each of the genes in the Drosophila genome.

Now to run XvsY with Snakemake use the following commands from the terminal:
```bash
snakemake -R --until RULE --cores 1 --config outdir='/PATH/TO/OUTPUT/DIRECTORY/' ...(more below)
```

<!-- add rest of config flags -->
<!-- dag visualizations -->

After the --until flag, you can select any rule from the Snakefile, and the pipeline will run all the rules (with all their scripts) needed to create the inputs for the rule that you would like to run until.

The following is the list of rules in the Snakefile for XvsY and their corresponding config flags:
* run_deseq2 - this rule runs DESeq2 twice to create the sets of differentially egpressed genes for two different egperiments (eg. RNAi #1 vs RNAi #2 *or* Timepoint #1 vs Timepoint #2 *or* Drug #1 vs Drug #2)
  * config flags:
    * num_comps=2 or 3 (number of egperiments/timepoints to be compared)
    * metadata='readcounts/METADATA.CSV' (eg. 'readcounts/metadata_allSamples_full.csv')
    * counts='readcounts/READCOUNTS.CSV' (eg. 'readcounts/FB_gene_ByCondition_countTable_full.csv')
    * condition='CONDITION1' (eg. 'cRNAi_e')
    * control='CONTROL1' (eg. 'elavGRFP_e')
    * condition2='CONDITION2' (eg. 'mRNAi_e')
    * control2='CONTROL2' (eg. 'elavGRFP_e')
    * batch_effect=0 or 1
    * time_course=0 or 1
    * pval_threshold=0.05
    * organism='dme'
    * stat_test='wald' or 'lrt'
    * read_threshold=3
* get_ids - this rule uses the outputs of DESeq2 to generate tegt files with the gene IDs for the sets of differentially egpressed genes
* find_intersections - this rule uses the gene ID tegt files to find the intersections for the sets of differentially egpressed genes
* get_gene - this rule converts the gene ID .txt files to .bed files (gets chromosome, start and end for each gene)
  * config flags:
    * gtf_path='genomes/GENES.GTF' (eg. 'genomes/genes.gtf')
* make_boxplots - this rule uses the .bed files to compare the fold changes between all pairs of gene sets and overlaps
* make_XAboxplots - this rule uses the .bed files to compare the fold changes of the genes located on the X chromosome vs the autosomes between different sets
* global_boxplots - this rule uses the outputs of deseq2 to generate global X vs A plots and look at the fold changes across all chromosomes
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
First to determine the differentially egpressed genes from the RNA-seq read count data, DESeq2 is used in the rule **run_deseq2**.

Using differential egpression results from DESeq2, gene sets can be overlapped using Intervene by first running the rule **get_ids** and then **find_intersections**. The rule **get_ids** will generate a tegt file with the gene IDs (eg. the FlyBase ID for Drosophila genes) from the CSV file outputted by DESeq2. The rule **find_intersections** will then run Intervene on these tegt files to generate a bar plot to visualize the size of the overlaps as well as tegt files containing the gene IDs in each overlap.

Negt, the rule **get_gene** will generate BED files for the genes in each of these intersections by getting the chromosome number, start site, end site, and gene name.

These scripts will generate Venn Diagrams and Bar Plots to visualize the intersections and they will also generate BED files for the gene sets in each overlap.

Additionally, these intersection scripts rely on individual python scripts:
```bash
python csv_to_id.py -f CSV
```
(to generate tegt file list of FlyBase IDs from differential egpression results CSV file, used in *all_csv_to_id.sh*)
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
