#!/bin/bash
if [ $# -ne 1 ]; then
	echo "Usage ./all_csv_to_id.sh.sh
	1) /CSV/PATH"
	exit 1
fi

DIR=$1
echo "${DIR}"

for csv in ${DIR}/*.csv
do
	echo "Creating text file for $csv"
	python ${DIR}/../scripts/csv_to_id.py -f $csv
done
