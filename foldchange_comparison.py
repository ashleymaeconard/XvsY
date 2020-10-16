import pandas as pd
import numpy as np
import argparse

import warnings
warnings.filterwarnings('ignore')

# create the command line parser
parser = argparse.ArgumentParser()

parser.add_argument("-m", "--file1", default='test1.csv', help="path to genes csv file to get data")
parser.add_argument("-n", "--file2", default='test2.csv', help="path to genes csv file to get data")


args = parser.parse_args()
data1 = args.file1
data2 = args.file2

# read data from csv file
df1 = pd.read_csv(data1)
df2 = pd.read_csv(data2)

# set genes with negative log2FoldChange to -1 and positive fold change to +1
df1[df1['log2FoldChange'] < 0] = -1
df1[df1['log2FoldChange'] > 0] = 1

df2[df2['log2FoldChange'] < 0] = -1
df2[df2['log2FoldChange'] > 0] = 1

df1 = df1.rename(columns = {'log2FoldChange':'log2FoldChange1'}, inplace = True)
df2 = df2.rename(columns = {'log2FoldChange':'log2FoldChange2'}, inplace = True)

# joins dataframes
df = df1.set_index('gene_id').join(df2.set_index('gene_id')).reset_index()

# keep rows not NA
df = df[df['log2FoldChange1'].notna()]
df = df[df['log2FoldChange2'].notna()]

# sum across rows
df['sum'] = df.sum(axis=1)

# count 0s
num_zeros = (df['sum'] != 0).sum()

print(num_zeros)