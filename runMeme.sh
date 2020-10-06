#!/bin/bash
if [ $# -ne 1 ]; then
	echo "Usage: ./runMeme.sh
	1) /PATH/TO/FOLDERS/WITH/GENELISTS (CLUSTERS)"
	exit 1
fi

#Inputting arguments
DIR=$1

#Checking all folders in directory for gene lists
for fol in "$DIR"/*
do
	echo "Checking directory: '$fol'"
	#Checking for CSVs in folder
	for csv in "$fol"/*
	do
		if [ ${csv: -4} == ".csv" ] && [ ! -d "$fol/MEME" ]; then
			echo "Checking GeneList: '$csv'"
			# Creating MEME input list of DNA sequences for set of genes in CSV
			python meme_suite_prep_indiv_cluster.py $DIR/../../reformatted_genes_gtf.csv $fol $DIR/../../dm6.fa 1 dme
			# Running MEME to identify motifs de novo for gene set and FIMO to find CLAMP motifs in gene set
			./run_meme_fimo_indiv_cluster.sh $fol $DIR/../../meme_CLAMP_overlap_GAF.txt
		fi
	done
done
