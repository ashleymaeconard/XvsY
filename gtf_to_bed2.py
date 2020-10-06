import pandas as pd
import numpy as np
import argparse
import os.path

import warnings
warnings.filterwarnings('ignore')

def create_gtf_bed(datapath):
    # load data from genes.gtf into dataframe
    df2 = pd.read_csv(datapath, sep="\t", names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    # select rows of dataframe for exons of genes
    exons = df2.loc[df2['feature'] == 'exon']
    # type conversions
    exons['attribute'] = exons['attribute'].astype(str)
    exons['seqname'] = exons['seqname'].astype(str)
    exons['start'] = exons['start'].astype(int)
    exons['end'] = exons['end'].astype(int)

    # getting data from the exons dataframe to create a genes bed file, which 
    # includes: seqname (chr), start, end, score, gene ids, and gene names
    genes = exons[['seqname', 'start', 'end']]
    exons['score'] = 500
    col5 = exons['score']
    split_atts = exons['attribute'].str.split("; ", expand = True) 
    geneIDs = split_atts[4].str.split(" ", n=1, expand = True) 
    geneIDs = geneIDs.loc[geneIDs[0] == 'gene_id']
    geneIDs = geneIDs[1].str.strip('\"')
    geneNames = split_atts[5].str.split(" ", n=1, expand = True) 
    geneNames = geneNames.loc[geneNames[0] == 'gene_name']
    geneNames = geneNames[1].str.strip('\"')

    # combine into new dataframe, finaldf
    finaldf = pd.DataFrame()
    finaldf = finaldf.join(genes, how='outer')
    finaldf['gene_id'] = geneIDs
    finaldf['score'] = col5
    finaldf['gene_name'] = geneNames

    # get a list of all unique genes in dataframe
    allgenes = finaldf['gene_name'].unique()

    # remove duplicate genes by combining exons into full length gene
    no_rep_df = pd.DataFrame(columns=['seqname', 'start', 'end', 'gene_id', 'score', 'gene_name'], index=None)
    for ngene in allgenes:
        dfgene = finaldf.loc[finaldf['gene_name'] == ngene]
        start = dfgene['start'].unique()
        start.sort()
        fst = start[0]
        end = dfgene['end'].unique()
        end = sorted(end, reverse=True)
        lst = end[0]
        dfgene['start'] = fst
        dfgene['end'] = lst
        dfgene.drop_duplicates(keep = 'last', inplace = True)
        no_rep_df = no_rep_df.append(dfgene)

    # rename dataframe columns for bed file
    no_rep_df.rename(columns = {'seqname':'chr'}, inplace = True)
    no_rep_df['chr'] = no_rep_df['chr'].astype(str)
    no_rep_df['start'] = no_rep_df['start'].astype(int)
    no_rep_df['end'] = no_rep_df['end'].astype(int)
    no_rep_df['gene_id'] = no_rep_df['gene_id'].astype(str)
    no_rep_df['score'] = no_rep_df['score'].astype(int)
    no_rep_df['gene_name'] = no_rep_df['gene_name'].astype(str)

    # save dataframe with no gene repeats, no_rep_df, as a bed file of all genes, gtf_genes.bed
    no_rep_df.to_csv(r'gtf_genes.bed', header=None, index=None, sep="\t", mode='w')
    print("Complete!")


def main():
    # create the command line parser
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--filepath", default='genes.gtf', help="path to genes.gtf file to get data")

    args = parser.parse_args()
    data = args.filepath

    #if gtf_genes.bed file already exists, do nothing, else run the create_gtf_bed function
    if os.path.exists("gtf_genes.bed"):
        print('GTF genes BED file already exists')
        return
    else:
        print('Creating BED file for GTF genes')
        create_gtf_bed(data)


if __name__ == '__main__':
    main() #runs main function
