#!/bin/bash
#Script to find intersections between TWO or THREE sets of genes and visualize them in venn diagrams and bar plots
#based on the FlyBase ID of these genes, and generates BED files each of these unique intersections.
#For up and downregulated gene sets COMBINED (not separate!)
if [ $# -ne 1 ]; then
	echo "Usage: ./intersections_combined_regulated_GeneSets.sh
	1) /PATH/TO/CSVs" #directory with all csv files (DESeq2 output matrix data)
	exit 1
fi

#Inputting arguments
DIR=$1

echo "Searching ${DIR}"

#Create txt files with FlyBase ID for each CSV experiment file
COUNTER=0
for csv in ${DIR}/*.csv
do
	echo "Creating txt file for $csv"
	python csv_to_id.py -f $csv
	COUNTER=$[COUNTER + 1]
done

#Create venn diagram
if [[ $COUNTER -eq 3 ]]
then
intervene venn -i ${DIR}/*.txt --type list --colors=green,darkorange,purple --figtype png --fontsize 12 --save-overlaps -o "$DIR/vennDiagram"
else
intervene venn -i ${DIR}/*.txt --type list --colors=green,purple --figtype png --fontsize 12 --save-overlaps -o "$DIR/vennDiagram"
fi

#Create bar plot
intervene upset -i ${DIR}/*.txt --type list --save-overlaps --figtype png -o "$DIR/barPlot"

#Create BED from GTF to get gene names, chromosome, start and end
python gtf_to_bed2.py -f ${DIR}/genes.gtf

#Get genes from FlyBase IDs in unique sets in the Venn Diagram and Bar Plot
for txv in ${DIR}/vennDiagram/sets/*.txt
do
	echo "In vennDiagram. Getting genes for $txv"
	python id_to_gene3.py -t $txv
done

for txb in ${DIR}/barPlot/sets/*.txt
do
	echo "In barPlot. Getting genes for $txb"
	python id_to_gene3.py -t $txb
done