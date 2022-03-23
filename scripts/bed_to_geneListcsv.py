# bed_to_geneListcsv.py
# Purpose: convert .bed to gene list
# Last mod. 02/22/2022

import pandas as pd
import numpy as np
import argparse
import glob
import os

import warnings
warnings.filterwarnings('ignore')

# create the command line parser
parser = argparse.ArgumentParser()

parser.add_argument("-f", "--filepath", default='010001_clampi_e_upreg_msl2i_e_upreg.bed', help="path to bed file to get gene list data")
parser.add_argument("-s", "--savepath", default='.', help="path to save gene list csv file")

args = parser.parse_args()
data = args.filepath
savepath = args.savepath

# read data from bed file into csv
df = pd.read_csv(data, sep="\t", names=["chr","start","end","gene_id","score","gene_name"])
rows,cols=df.shape

# get the fifth and fourth columns, gene_name and score respectively
df2 = df.iloc[:,[5,4]]
df2.rename(columns = {'gene_name':''}, inplace = True)
# set all scores to 1
df2["score"] = 1

# get filename to save df2
splitname = data.split("/")[-1].replace('.', '_').split("_")
filename = "_".join(splitname[:-1])
filename2 = filename + "_geneList_1.csv"

# get path to save df2
# cwd = os.getcwd()
# savepath = cwd+"/test_results/clusters/"+filename
savedir = savepath + "/" + filename

# create directory if it doesnt already exist
if not os.path.isdir(savedir):
	os.makedirs(savedir)
print("Directory '% s' created", savedir)

# save dataframe, df2
df2.to_csv(savedir+"/"+filename2, header=True, index=None, mode='w')
print("CSV Gene List Successfully Made! => ", filename2)
