import pandas as pd
import numpy as np
import argparse
import glob
import seaborn
from statannot import add_stat_annotation
import matplotlib.pyplot as plt
import itertools

import warnings
warnings.filterwarnings('ignore')

def makeViolinPlot(file_lst):
    li = []
    intersection = False
    for filename in file_lst:
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
            if 'e' in splitname:
                file1 = "_".join(splitname[2:5])
            else:
                file1 = "_".join(splitname[2:6])
            # read csv with expression data into dataframe and get gene fold changes
            mag_fold1 = pd.read_csv(file1 + ".csv")
            mag_fold1 = mag_fold1[['log2FoldChange']].astype(float)
            # add stuff into df
            df['log2FoldChange'] = mag_fold1
            df['dataset'] = file1
            df['overlap'] = file1 #+ ' (n:' + n + ')'
            # append dataframe to overall list for later visualization
            li.append(df)
        else:   # for intersection of genes between two experimental conditions, ie. overlapping group
            # read bed file into dataframe
            intersection = True
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
            df = pd.read_csv(filename, sep="\t", names=['chr', 'start', 'end', 'gene_id', 'score', 'gene_name'])
            n = str(len(df.index))
            mag_fold1 = pd.read_csv(file1 + ".csv")
            mag_fold1 = mag_fold1[['log2FoldChange']].astype(float)
            # add stuff to df
            df['log2FoldChange'] = mag_fold1
            df['dataset'] = file1
            df['overlap'] = fullfile #+ ' (n:' + n + ')'
            # append dataframe for first set to overall list for later visualization
            li.append(df)
            # read csv with expression data for first experimental condition into dataframe and get gene fold changes
            df = pd.read_csv(filename, sep="\t", names=['chr', 'start', 'end', 'gene_id', 'score', 'gene_name'])
            n = str(len(df.index))
            mag_fold2 = pd.read_csv(file2 + ".csv")
            mag_fold2 = mag_fold2[['log2FoldChange']].astype(float)
            # add stuff to df
            df['log2FoldChange'] = mag_fold2
            df['dataset'] = file2
            df['overlap'] = fullfile #+ ' (n:' + n + ')'
            # append dataframe for first set to overall list for later visualization
            li.append(df)

    # concatenate list into single larger dataframe
    frame = pd.concat(li, axis=0, ignore_index=True)

    # generate empty plot 
    plt.figure(figsize=(8,5))
    unique_datasets = frame['dataset'].unique()

    # determine up and down regulated gene sets
    done_up = False
    done_down = False
    for d in unique_datasets:
        if 'upreg' in d and not done_up:
            upreg_set = d
            done_up = True
        elif 'downreg' in d and not done_down:
            downreg_set = d
            done_down = False

    # my_pal = {"clampi_e_upreg (n:456)": "#af8dc3", "clampi_e_downreg (n:377)":"#7fbf7b"}
    # my_pal = {"msl2i_A_m_upreg": "#af8dc3", "msl2i_A_m_downreg": "#7fbf7b"}

    # set colors for different groups
    my_pal = {upreg_set: "#af8dc3", downreg_set: "#7fbf7b"}

    # get all unique datasets (experimental conditions) from dataframe, frame
    unique_datasets = frame['dataset'].unique()
    # create combinations of size 2
    subsets3 = list(itertools.combinations(unique_datasets, 2))
    box_pairs_lst = []
    for pair in subsets3:
        tup1,tup2 = pair
        ds1,chr1 = tup1
        ds2,chr2 = tup2
        if ds1 == ds2 or chr1 == chr2:
            box_pairs_lst.append(pair)

    # get all possible subsets of these combinations
    all_ss = []
    for L in range(1, len(box_pairs_lst)+1):
        for ss in itertools.combinations(box_pairs_lst, L):
            all_ss.append(list(ss))

    # sort subsets for largest to smallest by length
    all_ss_sorted = sorted(all_ss, key=len, reverse=True)

    # for each of these subsets try to create the statistical significance bars between different groups 
    # using the Mann-Whitney test
    for blabla in all_ss_sorted:
        try:
            if not intersection:
                # make violin plot of dataset (experimental condition) vs log2FoldChange
                plt.figure(figsize=(8,5))
                g = seaborn.violinplot(x="dataset", y="log2FoldChange", data=frame, palette=my_pal)
                # add significance bars
                add_stat_annotation(g, data=frame, x="dataset", y="log2FoldChange",
                                box_pairs=blabla,
                                test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
            elif intersection:
                # make violin plot of dataset (experimental condition) vs log2FoldChange
                g = seaborn.violinplot(x="dataset", y="log2FoldChange", data=frame, hue="overlap", palette=my_pal)
                plt.legend(title='Groups with Same Genes', bbox_to_anchor=(1.05,1), loc=2)
                # g.set_xticklabels(g.get_xticklabels(), rotation=20)
                # add significance bars
                add_stat_annotation(g, data=frame, x="dataset", y="log2FoldChange",
                                    box_pairs=blabla,
                                    test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
            break
        except ValueError:
            continue

    plt.tight_layout()

    plt.tick_params(axis='x', pad=17)

    # add sample size (n: _) under each of the violin groups
    medians = frame.groupby(['dataset'])['log2FoldChange'].median() 
    mins = frame.groupby(['dataset'])['log2FoldChange'].min() 
    nobs =  frame.groupby(['dataset']).apply(lambda x: 'n: {}'.format(len(x)))
    ymin,ymax = g.get_ylim()
    for ax in plt.gcf().axes:
        for tick, label in enumerate(ax.get_xticklabels()):
            ax_dat = label.get_text()
            x_offset = 0
            min_val = mins[ax_dat]
            num = nobs[ax_dat]
            ax.text(tick + x_offset, ymin - abs(0.0686*(ymax-ymin)), num,  #ymin - 0.25
                    horizontalalignment='center', size='medium', color='green', weight='semibold')

    # get title and save figure
    fig = g.get_figure()
    if 'e' in splitname:
        figtit = "_".join(splitname[2:5])
    else:
        figtit = "_".join(splitname[2:6])
    fig.savefig('overall_violin_'+figtit+'.png')



def main():
    # create the command line parser
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--filepath", default='hello.bed', help="path to bed file to get data")

    args = parser.parse_args()
    data = args.filepath

    # get list of bed files with genes sets to compare
    all_files = glob.glob(data + "/*.bed")

    for idx,f in enumerate(all_files):
        makeViolinPlot(f)
        remaining_files = all_files[idx:]
        for f2 in remaining_files:
            file_lst = [f,f2]
            makeViolinPlot(file_lst)


if __name__ == "__main__":
    main()