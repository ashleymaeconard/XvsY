#!/bin/bash
if [ $# -ne 2 ]; then
	echo "Usage: ./runGOall.sh
		1) /PATH/TO/BED_DIR
		2) /PATH/TO/INPUT_OUTPUT_DIR/ "
	exit 1
fi

#Inputting arguments
BED_DIR=$1
OUT_DIR=$2

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
	Rscript scripts/clusterProfiler.r $godir $actualdir 0 dme 0.05
done
echo "Done!"
