import pandas as pd
import numpy as np
import argparse
import os

import warnings
warnings.filterwarnings('ignore')

# create the command line parser
parser = argparse.ArgumentParser()

parser.add_argument("-f", "--filepath", default='clusters', help="path to get data")
parser.add_argument("-g", "--gogroup", default='BP', help="BP, MF, CC")

args = parser.parse_args()
rootdir = args.filepath
goterm = args.gogroup

data_df = pd.DataFrame(columns=['Set', 'Description'], index=None)

if rootdir == '.':
    rootdir = os.getcwd()

if goterm == 'BP':
    others = ['MF','CC']
if goterm == 'MF':
    others = ['BP','CC']
if goterm == 'CC':
    others = ['BP','MF']

for root, subdirs, files in os.walk(rootdir):
    for subdir in subdirs:
        filename = subdir
        for root1, subdirs1, files1 in os.walk(rootdir+'/'+subdir):
            if subdirs1 == []:
                for f in files1:
                    if '.csv' in f:
                        set_name = filename.split('_')
                        set_name = '_'.join(set_name[1:])
                        int_df = pd.DataFrame(columns=['Set', 'Description'], index=None)
                        int_df.loc[0, 'Set'] = set_name
                        int_df.loc[0, 'Description'] = 'EMPTY'
                        data_df = data_df.append(int_df)
            if goterm in subdirs1:
                subdirs1 = [goterm]
            multiples = True
            for subdir1 in subdirs1:
                if goterm in subdir1:
                    df = pd.read_csv(rootdir+'/'+subdir+'/'+subdir1+'/'+filename+'_cluster__1_clustProf_'+goterm+'_.csv')
                    final_df = df.sort_values('pvalue')
                    vals = final_df.iloc[0, [2,5]]
                    top_go_term = final_df.iloc[0,2]
                    set_name = filename.split('_')
                    set_name = '_'.join(set_name[1:])
                    int_df = pd.DataFrame(columns=['Set', 'Description'], index=None)
                    int_df.loc[0, 'Set'] = set_name
                    int_df.loc[0, 'Description'] = top_go_term
                    data_df = data_df.append(int_df)
                elif subdir1 in others and multiples:
                    set_name = filename.split('_')
                    set_name = '_'.join(set_name[1:])
                    int_df = pd.DataFrame(columns=['Set', 'Description'], index=None)
                    int_df.loc[0, 'Set'] = set_name
                    int_df.loc[0, 'Description'] = 'EMPTY'
                    data_df = data_df.append(int_df)
                    multiples = False
                elif subdir1 == 'www':
                    pass
                    
data_df.drop_duplicates(keep = 'last', inplace = True)
print(data_df)
if not os.path.isdir(rootdir+'/../summary/'):
    os.makedirs(rootdir+'/../summary/')
data_df.to_csv(rootdir+'/../summary/'+goterm+'_set_terms.csv',index=None)
