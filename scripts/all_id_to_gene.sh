#!/bin/bash
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

python ${DIR}/../../../scripts/gtf_to_bed2.py -f ${GENESGTF} -s $GENESBED

for txf in ${DIR}/*.txt
do
	echo "Creating BED file for $txf"
	if [[ "${DIR}" == *"$SUB"* ]]; 
	then
		python ${DIR}/../../../scripts/id_to_gene3.py -t $txf -b $GENESBED
	else
		python ${DIR}/../../scripts/id_to_gene3.py -t $txf -b $GENESBED
	fi
	
done
