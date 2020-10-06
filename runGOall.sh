#!/bin/bash
if [ $# -ne 1 ]; then
	echo "Usage: ./runGOall.sh
		1) /PATH/TO/BED_DIR"
	exit 1
fi

#Inputting arguments
DIR=$1

echo "Searching ${DIR}"

for FILE in ${DIR}/*.bed
do
	echo "Getting Gene List CSV"
	python bed_to_geneListcsv.py -f $FILE

	#Getting DIR name
	savedir="${FILE%.*}"
	IFS='/' read -r -a array <<< $savedir
	actualdir="${array[-1]}"
	#echo $actualdir
	BASEDIR="$(dirname $FILE)"
	#echo $BASEDIR
	med="/test_results/clusters/"
	godir="$BASEDIR$med$actualdir/"
	#echo $godir

	echo "Running GO Analysis the Gene List"
	Rscript clusterProfiler.r $godir $actualdir 0 dme 0.05
done
echo "Done!"
