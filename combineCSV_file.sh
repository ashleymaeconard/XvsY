#!/bin/bash
#Script to combine two specific CSV files into one with a given name
if [ $# -ne 3 ]; then
	echo "Usage: ./combineCSV.sh
		1) /PATH/TO/CSV/1
		2) /PATH/TO/CSV/2
		3) Output File Name"
	exit 1
fi

#Inputting arguments
csv1=$1
csv2=$2
outputfile=$3

#Combining CSVs
head -1 "$csv1" > "$outputfile"
tail -n +2 "$csv1" >> "$outputfile"
tail -n +2 "$csv2" >> "$outputfile"
