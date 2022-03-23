#!/bin/bash
# runGOall.sh
# Purpose: prepare inputs to run GO analysis with clusterProfiler
# Last mod. 02/22/2022

if [ $# -ne 5 ]; then
	echo "Usage: ./runGOall.sh
		1) /PATH/TO/BED_DIR
		2) /PATH/TO/INPUT_OUTPUT_DIR/ 
		3) SEPARATE_TIMEPOINTS (set to 1 to run GO for each timepoint separately, 0 otherwise; 0 recommended as default)
		4) ORGANISM (dme, hsa, mmu)
		5) ADJUSTED_PVALUE (recommend 0.05)"
	exit 1
fi

#Inputting arguments
BED_DIR=$1
OUT_DIR=$2
SEP_TPS=$3
ORGANISM=$4
ADJ_PVAL=$5

echo "Searching ${BED_DIR}"

for FILE in ${BED_DIR}/*.bed
do
	echo "Getting Gene List CSV"
	echo ${BED_DIR}
	echo $FILE
	python scripts/bed_to_geneListcsv.py -f $FILE -s $OUT_DIR

	#Getting DIR name
	savedir="${FILE%.*}"
	IFS='/' read -r -a array <<< $savedir
	actualdir="${array[-1]}"
	echo $actualdir
	# BASEDIR="$(dirname $FILE)"
	#echo $BASEDIR
	# med="/test_results/clusters/"
	# godir="$BASEDIR$med$actualdir/"
	#echo $godir
	slash="/"
	godir="$OUT_DIR$slash$actualdir"
	echo $godir

	echo "Running GO Analysis the Gene List"
	Rscript scripts/clusterProfiler.r $godir $actualdir $SEP_TPS $ORGANISM $ADJ_PVAL
done
echo "Done!"
