#!/bin/bash
if [ $# -ne 6 ]; then
	echo "Usage: ./runMEMEall.sh
		1) /PATH/TO/BED_DIR
		2) /PATH/TO/INPUT_OUTPUT_DIR/
        3) /PATH/TO/reformatted_genes_gtf.csv
        4) /PATH/TO/CHROM_FAs
        5) TSS_only (set to 1 to run +-1kb from transcription start site, 0 otherwise)
        6) GENOME (dm6, hsa, or mmu)"
	exit 1
fi

#Inputting arguments
BED_DIR=$1
OUT_DIR=$2
REFORM_GENES=$3
CHROM_FA=$4
TSS_ONLY=$5
ORGANISM=$6

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
                python scripts/meme_suite_prep_indiv_cluster.py $REFORM_GENES $fol $CHROM_FA $TSS_ONLY $ORGANISM
                echo "FASTA made!"
            fi
        done
    done

done
echo "Done!"
