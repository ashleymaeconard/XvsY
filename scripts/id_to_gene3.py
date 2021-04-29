import pandas as pd
import numpy as np
import argparse
#import sys
import os

import warnings
warnings.filterwarnings('ignore')

# create the command line parser
parser = argparse.ArgumentParser()

parser.add_argument("-t", "--txtfile", default='test.txt', help="path to id .txt file to get data")
parser.add_argument("-b", "--bedfile", default='gtf_genes.bed', help="path to genes gtf .bed file")

args = parser.parse_args()
ids = args.txtfile
genes_pth = args.bedfile

if (genes_pth[-1] != '/'):
	genes_pth = genes_pth + '/'

#ids = sys.argv[1]

if not os.path.isdir(os.getcwd()+'/deg_sets/gene_beds'):
	os.mkdir(os.getcwd()+'/deg_sets/gene_beds')

# opens text files of FlyBase IDs that are in a unique intersection
flybase_ids_txt = open(ids, 'r')
# create list of all FlaseBase IDs and strip new line character
lst_ids = [line.strip('\n') for line in flybase_ids_txt.readlines()]

# read gtf_genes.bed file as dataframe
genes = pd.read_csv(genes_pth+"gtf_genes.bed", sep="\t", names=['chr', 'start', 'end', 'gene_id', 'score', 'gene_name'])

# find genes present in both the gtf_genes.bed file as well as the gene list from the intersection text file
overlaps = genes[genes['gene_id'].isin(lst_ids)]

# create new dataframe with the chr, start, end, gene_id, score, and gene_name of these common genes
bed = pd.DataFrame()
bed['chrom'] = overlaps['chr'].astype(str)
bed['chromStart'] = overlaps['start'].astype(int)
bed['chromEnd'] = overlaps['end'].astype(int)
bed['geneID'] = overlaps['gene_id'].astype(str)
bed['score'] = overlaps['score'].astype(int)
bed['geneName'] = overlaps['gene_name'].astype(str)

dros_chroms = ['X', '2R', '2L', '3R', '3L', '4', 'Y']
bed = bed.loc[bed['chrom'].isin(dros_chroms)]

# save dataframe as bed file
name = ids.split('/')[-1].split(".")[0] + ".bed"
bed.to_csv(os.getcwd()+'/deg_sets/gene_beds/'+name, header=None, index=None, sep="\t", mode='w')
print('Done!')
