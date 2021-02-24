#!/bin/bash
if [ $# -ne 2 ]; then
	echo "Usage: ./runMEMEall.sh
		1) /PATH/TO/BED_DIR
		2) /PATH/TO/INPUT_OUTPUT_DIR/"
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

	for fol in "$OUT_DIR"/*
    do
        #echo "Checking directory: '$fol'"
        for csv in "$fol"/*
        do
            if [ ${csv: -4} == ".csv" ] && [ ! -d "$fol/MEME" ]; then
                echo "Checking GeneList: '$csv'"
                python scripts/meme_suite_prep_indiv_cluster.py $OUT_DIR/../../scripts/reformatted_genes_gtf.csv $fol $OUT_DIR/../../scripts/dm6.fa 1 dme
                echo "FASTA made!"
            fi
        done
    done

done
echo "Done!"
