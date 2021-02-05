import pandas as pd
import numpy as np
import argparse
import glob
import seaborn
from statannot import add_stat_annotation
import matplotlib.pyplot as plt
import itertools
import 

import warnings
warnings.filterwarnings('ignore')

def makeViolinPlot(genegroup_lst, savepath):

    li = []

    for filename in genegroup_lst:
        temp = filename.replace('.', '_').split("_") #[2:-1]
        print(temp)
        # for genes unique to a specific experimental condition, ie. non-overlapping group
        if len(temp) < 8:
            # read bed file into dataframe
            df = pd.read_csv(filename, sep="\t", names=['chr', 'start', 'end', 'gene_id', 'score', 'gene_name'])
            n = str(len(df.index))
            splitname = filename.replace('.', '_').split("_")
            print(splitname)
            # get file base name, shorter if embryo which is unsexed
            if 'e' in temp:
                fullfile = "_".join(splitname[2:5])
                file1 = "_".join(splitname[2:5])
            else:
                fullfile = "_".join(splitname[2:6])
                file1 = "_".join(splitname[2:6])
            # read csv with expression data into dataframe and get gene fold changes
            mag_fold1 = pd.read_csv(file1 + ".csv")
            gene_ids = df['gene_id'].tolist()
            mag_fold1 = mag_fold1[mag_fold1['log2FoldChange'].isin(gene_ids)].astype(float)
            # add stuff into df
            df['log2FoldChange'] = mag_fold1
            df['dataset'] = file1
            df['overlap'] = fullfile
            # change all chromosome labels to either X or A (autosome)
            df['chr'] = df['chr'].map({'2L':'A', '2R':'A', '3L':'A', '3R':'A', '4':'A', 'X':'X'})
            # append dataframe to overall list for later visualization
            li.append(df)
        else:   # for intersection of genes between two experimental conditions, ie. overlapping group
            # read bed file into dataframe
            df = pd.read_csv(filename, sep="\t", names=['chr', 'start', 'end', 'gene_id', 'score', 'gene_name'])
            n = str(len(df.index))
            splitname = filename.replace('.', '_').split("_")
            print(splitname)
            # get base names of each files, shorter if embryo which is unsexed
            if 'e' in temp:
                fullfile = "_".join(splitname[2:8])
                file1 = "_".join(splitname[2:5])
                file2 = "_".join(splitname[5:8])
            else:
                fullfile = "_".join(splitname[2:10])
                file1 = "_".join(splitname[2:6])
                file2 = "_".join(splitname[6:10])
            # read csv with expression data for first experimental condition into dataframe and get gene fold changes
            mag_fold1 = pd.read_csv(file1 + ".csv")
            mag_fold1 = mag_fold1[['log2FoldChange']].astype(float)
            # add stuff to df
            df['log2FoldChange'] = mag_fold1
            df['dataset'] = file1
            df['overlap'] = fullfile
            # change all chromosome labels to either X or A (autosome)
            df['chr'] = df['chr'].map({'2L':'A', '2R':'A', '3L':'A', '3R':'A', '4':'A', 'X':'X'})
            # append dataframe for first set to overall list for later visualization
            li.append(df)
            # read csv with expression data for second experimental condition into dataframe and get gene fold changes
            df = pd.read_csv(filename, sep="\t", names=['chr', 'start', 'end', 'gene_id', 'score', 'gene_name'])
            n = str(len(df.index))
            mag_fold2 = pd.read_csv(file2 + ".csv")
            mag_fold2 = mag_fold2[['log2FoldChange']].astype(float)
            # add stuff to df
            df['log2FoldChange'] = mag_fold2
            df['dataset'] = file2
            df['overlap'] = fullfile
            # change all chromosome labels to either X or A (autosome)
            df['chr'] = df['chr'].map({'2L':'A', '2R':'A', '3L':'A', '3R':'A', '4':'A', 'X':'X'})
            # append dataframe for second set to overall list for later visualization
            li.append(df)

    # concatenate list into single larger dataframe
    frame = pd.concat(li, axis=0, ignore_index=True)

    # get all unique datasets (experimental conditions) from dataframe, frame
    unique_datasets = frame['dataset'].unique()
    xa_chrom = ['X', 'A']
    # make into tuples, ie ('X', clampi_e)
    subsets2 = list(itertools.product(unique_datasets, xa_chrom))
    print(subsets2)
    # create combinations of size 2 between pairs of tuples
    subsets3 = list(itertools.combinations(subsets2, 2))
    print(subsets3)
    box_pairs_lst = []
    for pair in subsets3:
        tup1,tup2 = pair
        ds1,chr1 = tup1
        ds2,chr2 = tup2
        if ds1 == ds2 or chr1 == chr2:
            box_pairs_lst.append(pair)
    print(box_pairs_lst)

    # get all possible subsets of these combinations
    all_ss = []
    for L in range(1, len(box_pairs_lst)+1):
        for ss in itertools.combinations(box_pairs_lst, L):
            all_ss.append(list(ss))
            print(ss)

    print(all_ss)

    # make violin plot of dataset (experimental condition) vs log2FoldChange with X vs A separated
    plt.figure(figsize=(8,5))
    g = seaborn.violinplot(x="dataset", y="log2FoldChange", data=frame, hue="chr")

    # plt.show()

    # sort subsets for largest to smallest by length
    all_ss_sorted = sorted(all_ss, key=len, reverse=True)
    print(all_ss_sorted[0])
    print(box_pairs_lst)

    # for each of these subsets try to create the statistical significance bars between different groups 
    # using the Mann-Whitney test
    for blabla in all_ss_sorted:
        print(blabla)
        try:
            plt.figure(figsize=(8,5))
            g = seaborn.violinplot(x="dataset", y="log2FoldChange", data=frame, hue="chr")
            add_stat_annotation(g, data=frame, x="dataset", y="log2FoldChange", hue="chr",
                            box_pairs=blabla,
                            test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
            break
        except ValueError:
            continue


    plt.tick_params(axis='x', pad=17)

    plt.tight_layout()

    # add sample size (n: _) under each of the violin groups
    medians = frame.groupby(['dataset','chr'])['log2FoldChange'].median() 
    mins = frame.groupby(['dataset','chr'])['log2FoldChange'].min() 
    nobs =  frame.groupby(['dataset','chr']).apply(lambda x: 'n: {}'.format(len(x)))
    ymin,ymax = g.get_ylim()
    for ax in plt.gcf().axes:
        for tick, label in enumerate(ax.get_xticklabels()):
            ax_dat = label.get_text()
            for j, ax_chr in enumerate(ax.get_legend_handles_labels()[1]):
                if len(ax.get_legend_handles_labels()[1]) == 1:
                    x_offset = 0
                else:
                    x_offset = (j - 0.5)*0.4
                min_val = mins[ax_dat, ax_chr]
                num = nobs[ax_dat, ax_chr]
                ax.text(tick + x_offset,  ymin - abs(0.0686*(ymax-ymin)), num,
                        horizontalalignment='center', size='medium', color='green', weight='semibold')

    fig = g.get_figure()

    # figtit = max(frame['overlap'].unique(), key=len)

    figtit = "_".join(frame['dataset'].unique())

    # save figure
    fig.savefig(savepath + '/XAviolin_'+figtit+'.png')

def main():
    # create the command line parser
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--filepath", default='.', help="path to bed file folder to get data")
    parser.add_argument("-s", "--savepath", default='.', help="path to save figures")

    args = parser.parse_args()
    data = args.filepath
    savepath = args.savepath

    # get list of bed files
    all_files = glob.glob(data + "/*.bed")

    if not os.path.isdir(os.getcwd()+'/comparisons/xa_boxplots'):
        os.mkdir(os.getcwd()+'/comparisons/xa_boxplots')

    for idx,f in enumerate(all_files):
        makeViolinPlot([f], savepath)
        remaining_files = all_files[idx+1:]
        for f2 in remaining_files:
            file_lst = [f,f2]
            makeViolinPlot(file_lst, savepath)


if __name__ == "__main__":
    main()