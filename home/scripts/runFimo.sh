#!/bin/bash
# runFIMO.sh
# Purpose: prepare inputs to run FIMO
# Reference: https://meme-suite.org/meme/doc/fimo.html
# Last mod. 02/22/2022

if [ $# -ne 2 ]; then
	echo "Usage: ./runMeme.sh
	1) /PATH/TO/FOLDERS/WITH/GENELISTS (CLUSTERS)
    2) /PATH/TO/PWM"
	exit 1
fi

#Inputting arguments
INPUT_DIR=$1
PWM=$2

for fol in "$INPUT_DIR"/*
do
	#echo "Checking directory: '$fol'"
	for geneList in "$fol"/MEME/*DNAseqs*
	do
		echo $fol
		echo $geneList
        echo $PWM
		# Setting output directory
    	OUTPUT_DIR=$(dirname "${geneList}")

        echo "Running FIMO"
    	fimo --oc $OUTPUT_DIR $PWM $geneList
	done
done
