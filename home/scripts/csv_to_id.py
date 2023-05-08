# csv_to_id.py
# Purpose: converts .csv to ids
# Last mod. 02/22/2022

import pandas as pd
import numpy as np
import argparse
import os
#import sys

import warnings
warnings.filterwarnings('ignore')

# create the command line parser
parser = argparse.ArgumentParser()

parser.add_argument("-f", "--filename", default='test.csv', help="path to genes csv file to get data")

args = parser.parse_args()
data = args.filename

#data = sys.argv[1]

print(data)

if not os.path.isdir(os.getcwd()+'/deg_sets/gene_ids'):
	os.mkdir(os.getcwd()+'/deg_sets/gene_ids')

# read data from csv file
df = pd.read_csv(data)

# filters genes by fold change p-values below 0.01
df = df.loc[df['padj'] < 0.01]

# checks if dataframe has gene_id column, if not just use first column
try: 
	ids = df[['gene_id']]
except:
	ids = df.iloc[:,0]

# saves FlyBase IDs for genes in csv as a text file
name = data.split(".")[0] + ".txt"
name = name.split("/")[-1]
print(name)
print(os.getcwd())
ids.to_csv(os.getcwd()+'/deg_sets/gene_ids/'+name, header=None, index=None, mode='w')
print("Success!")
