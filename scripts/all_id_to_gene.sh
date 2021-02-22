#!/bin/bash
if [ $# -ne 1 ]; then
	echo "Usage ./all_csv_to_id.sh.sh
	1) /GENE/TXT/PATH"
	exit 1
fi

DIR=$1
echo "${DIR}"
SUB='barPlot/set'

for txf in ${DIR}/*.txt
do
	echo "Creating BED file for $txf"
	if [[ "${DIR}" == *"$SUB"* ]]; 
	then
		python ${DIR}/../../../scripts/id_to_gene3.py -t $txf
	else
		python ${DIR}/../../scripts/id_to_gene3.py -t $txf
	fi
	
done
