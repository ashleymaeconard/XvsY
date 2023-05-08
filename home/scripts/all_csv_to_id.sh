#!/bin/bash
# all_csv_to_id.sh
# Purpose: call funciton to convert csvs to ids for all .csvs in input directory
# Last mod. 02/22/2022

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
