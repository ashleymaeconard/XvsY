## RNA-seq Analysis Pipeline for data across sex and time
CLAMP-null vs CLAMP RNAi vs MSL2 RNAi Effects on Brain Development
Tentative Name: EXPressST

## Intersections
Using differential expression results from DESeq2, gene sets can be overlapped using Intervene by running:
```bash
./intersections_separated_regulated_GeneSets.sh /PATH/TO/CSVs
```

The above command will work if the gene sets have been separated into upregulated and downregulated. To combine these gene sets into an overall group of genes that are differentially expressed in a given experimental condition, run:
```bash
./combineCSV_dir.sh /PATH/TO/CSVs
```
(to combine all CSVs in directory)
OR
```bash
./combineCSV_file.sh /PATH/TO/CSV1 /PATH/TO/CSV2 Output_File_Name
```
(to combine two specific files with a given output name)

Then the combined gene sets can be overlapped to find the intersections by running:
```bash
./intersections_combined_regulated_GeneSets.sh /PATH/TO/CSVs
```

These scripts will generate Venn Diagrams and Bar Plots to visualize the intersections and they will also generate BED files for the gene sets in each overlap.

Additionally, these intersection scripts rely on individual python scripts:
```bash
python csv_to_id.py -f CSV
```
(to generate text file list of FlyBase IDs from differential expression results CSV file)
```bash
python gtf_to_bed2.py -f /PATH/TO/genes.gtf 
```
(to create BED file from GTF to get gene names, chromosome, start and end)
```bash
python id_to_gene3.py -t TXT_FILE 
```
(to create BED files to Get genes from FlyBase IDs in unique sets from Intervene)

## Fold Change Comparison between gene sets
To compare the fold changes between two or more of these gene sets in BED file format (from the intersections), run:
```bash
python compareGeneSet_boxplots.py -f /PATH/TO/BED/FILES 
```

## X vs A Plots
To generate global X vs A boxplots for a given gene set in BED file format (from the intersections), run:
```bash
python global_XA.py -f BED_FILE -c CASE_ID
```

To generate violin plots comparing the X vs A fold change between two gene groups in BED file format, run:
```bash
python xa_boxplots.py -f /PATH/TO/BED/FILES 
```

# GO Analysis
To run gene ontology (GO) analysis to determine the enriched biological processes (BP), molecular function (MF), and cellular component (CC) of given gene set as a BED file, run:
```bash
./runGO.sh /PATH/TO/BED
```
OR to run GO analysis for all BED files in directory, use:
```bash
./runGOall.sh /PATH/TO/BED_DIR
```

These scripts use: bed_to_geneListcsv.py and clusterProfiler.r

## MEME and FIMO Analysis
To run MEME and FIMO analysis to find the top three de novo motifs in a specific gene set and the presence of CLAMP motifs respectively, run:
```bash
./runMeme.sh /PATH/TO/FOLDERS/WITH/GENELISTS_CSV (CLUSTERS)
```

This script uses: meme_suite_prep_indiv_cluster.py and run_meme_fimo_indiv_cluster.sh

To generate a summary of the FIMO results as a table showing the percentage of a each motif, run:
```bash
python meme_summary.py -p /PATH/TO/MEME_RESULTS
```

![illustration of pipline](https://github.com/ashleymaeconard/clampvsmsl2/blob/develop/workflow_diagram2.PNG?raw=true)