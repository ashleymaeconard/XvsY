#!/bin/bash
#Script to combine all CSV files from same experiment (ex. up and downregulated genes with RNAi) into one file
if [ $# -ne 1 ]; then
	echo "Usage: ./combineCSV.sh
		1) /PATH/TO/CSVs"
	exit 1
fi

#Inputting arguments
DIR=$1

#Combining CSVs
for csv1 in ${DIR}/*.csv
do
	file1="${csv1%_*}" #Gets first CSV file (ex. clampi_e_upreg.csv)
	for csv2 in ${DIR}/*.csv
	do
		file2="${csv2%_*}" #Gets second CSV file (ex. clampi_e_downreg.csv)
		if [ "$file1" == "$file2" -a "$csv1" != "$csv2" ] ;
		then
			name="${file1##*/}"
			outputfile="combined_$name.csv" #Name of new combined CSV file
			if [ ! -f "$outputfile" ] ;
			then
				echo $outputfile
				#Combines files with the header of one and 
				#everything but the first line of both
				head -1 "$csv1" > "$outputfile"
				tail -n +2 "$csv1" >> "$outputfile"
				tail -n +2 "$csv2" >> "$outputfile"
				#Removes old CSV files
				rm -f "$csv1"
				rm -f "$csv2"	
			fi
		fi
	done
done
