import pandas as pd
import numpy as np
import argparse
import os.path

import warnings
warnings.filterwarnings('ignore')

def create_gtf_bed(datapath,savepath):
    df2 = pd.read_csv(datapath, sep="\t", names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    exons = df2.loc[df2['feature'] == 'exon']
    exons['attribute'] = exons['attribute'].astype(str)
    exons['seqname'] = exons['seqname'].astype(str)
    exons['start'] = exons['start'].astype(int)
    exons['end'] = exons['end'].astype(int)

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

    finaldf = pd.DataFrame()
    finaldf = finaldf.join(genes, how='outer')
    finaldf['gene_id'] = geneIDs
    finaldf['score'] = col5
    finaldf['gene_name'] = geneNames
    allgenes = finaldf['gene_name'].unique()

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

    no_rep_df.rename(columns = {'seqname':'chr'}, inplace = True)
    no_rep_df['chr'] = no_rep_df['chr'].astype(str)
    no_rep_df['start'] = no_rep_df['start'].astype(int)
    no_rep_df['end'] = no_rep_df['end'].astype(int)
    no_rep_df['gene_id'] = no_rep_df['gene_id'].astype(str)
    no_rep_df['score'] = no_rep_df['score'].astype(int)
    no_rep_df['gene_name'] = no_rep_df['gene_name'].astype(str)

    no_rep_df.to_csv(savepath+'gtf_genes.bed', header=None, index=None, sep="\t", mode='w')
    print("Complete!")


def main():
    # create the command line parser
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--filepath", default='genes.gtf', help="path to genes.gtf file to get data")
    parser.add_argument("-s", "--savepath", default='.', help="path to save .bed file made from genes.gtf file")

    args = parser.parse_args()
    data = args.filepath
    savepath = args.savepath

    if (savepath[-1] != '/'):
        savepath = savepath + '/'

    if os.path.exists(savepath+"gtf_genes.bed"):
        print('GTF genes BED file already exists')
        return
    else:
        print('Creating BED file for GTF genes')
        create_gtf_bed(data,savepath)


if __name__ == '__main__':
    main()
