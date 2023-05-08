#!/bin/bash
# all_id_to_gene.sh
# Purpose: call funciton to convert ids to associated gene name (.bed format) for all .txts in input directory
# Last mod. 02/22/2022

if [ $# -ne 2 ]; then
	echo "Usage ./all_csv_to_id.sh.sh
	1) /GENE/TXT/DIR
	2) /GENES/GTF/PATH"
	exit 1
fi

DIR=$1
GENESGTF=$2
echo "${DIR}"
SUB='barPlot/set'
GENESBED=$(dirname "$GENESGTF")

python ${DIR}/../../../scripts/gtf_to_bed.py -f ${GENESGTF} -s $GENESBED

for txf in ${DIR}/*.txt
do
	echo "Creating BED file for $txf"
	if [[ "${DIR}" == *"$SUB"* ]]; 
	then
		python ${DIR}/../../../scripts/id_to_gene_map.py -t $txf -b $GENESBED
	else
		python ${DIR}/../../scripts/id_to_gene_map.py -t $txf -b $GENESBED
	fi
	
done
