#!/bin/bash
# runMeme.sh
# Purpose: prepare inputs to run MEME
# Reference: https://meme-suite.org/meme/
# Last mod. 02/22/2022

if [ $# -ne 1 ]; then
	echo "Usage: ./runMeme.sh
	1) /PATH/TO/FOLDERS/WITH/GENELISTS (CLUSTERS)"
	exit 1
fi

#Inputting arguments
INPUT_DIR=$1

for fol in "$INPUT_DIR"/*
do
	#echo "Checking directory: '$fol'"
	for geneList in "$fol"/MEME/*DNAseqs*
	do
		echo "Running MEME"
		echo $fol
		echo $geneList
		# Setting output directory
    	OUTPUT_DIR=$(dirname "${geneList}")
    	# -objfun classic --revcomp
    	meme $geneList -dna -nmotifs 1 -maxsize 2020000 -mod anr -oc $OUTPUT_DIR/	
	done
done
