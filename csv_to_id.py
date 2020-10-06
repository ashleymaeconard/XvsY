import pandas as pd
import numpy as np
import argparse

import warnings
warnings.filterwarnings('ignore')

# create the command line parser
parser = argparse.ArgumentParser()

parser.add_argument("-f", "--filename", default='test.csv', help="path to genes csv file to get data")

args = parser.parse_args()
data = args.filename

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
ids.to_csv(name, header=None, index=None, mode='w')
print("Success!")
