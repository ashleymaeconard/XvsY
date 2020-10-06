#!/bin/bash
if [ $# -ne 1 ]; then
	echo "Usage: ./runGO.sh
		1) /PATH/TO/BED"
	exit 1
fi

#Inputting arguments
FILE=$1

#Getting Gene List CSV
echo "Getting Gene List CSV"
python bed_to_geneListcsv.py -f $FILE

#Getting DIR name
savedir="${FILE%.*}"
IFS='/' read -r -a array <<< $savedir
actualdir="${array[-1]}"
echo $actualdir
BASEDIR="$(dirname $FILE)"
echo $BASEDIR
med="/test_results/clusters/"
godir="$BASEDIR$med$actualdir/"
echo $godir

#Running clusterProfiler for GO Analysis
echo "Running GO Analysis the Gene List"
Rscript clusterProfiler.r $godir $actualdir 0 dme 0.05

echo "Done!"
