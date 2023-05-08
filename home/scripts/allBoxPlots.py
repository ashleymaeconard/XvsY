# allBoxPlots.py
# Purpose: create boxplots
# Last mod. 02/22/2022

import pandas as pd
import numpy as np
import argparse
import glob
import seaborn
from statannot import add_stat_annotation
import matplotlib.pyplot as plt
import itertools
import os

import warnings
warnings.filterwarnings('ignore')

def makeViolinPlot(file_lst, degpath, savepath):
    li = []
    intersection = False
    for filename in file_lst:
        count_downreg1 = file_lst[0].count('downreg')
        count_upreg1 = file_lst[0].count('upreg')
        total_counts1 = count_downreg1 + count_upreg1
        print(total_counts1)
        count_downreg2 = file_lst[1].count('downreg')
        count_upreg2 = file_lst[1].count('upreg')
        total_counts2 = count_downreg2 + count_upreg2
        print(total_counts2)
        if total_counts1 > 2 or total_counts2 > 2:
            return
        temp = filename.split('/')[-1].replace('.', '_').split("_") #[2:-1]
        print(temp)
        # for genes unique to a specific experimental condition, ie. non-overlapping group
        if len(temp) < 8:
            # read bed file into dataframe
            df = pd.read_csv(filename, sep="\t", names=['chr', 'start', 'end', 'gene_id', 'score', 'gene_name'])
            n = str(len(df.index))
            splitname = filename.split('/')[-1].replace('.', '_').split("_")
            print(splitname)
            # get file base name, shorter if embryo which is unsexed
            if 'e' in splitname:
                file1 = "_".join(splitname[1:-1])
            else:
                file1 = "_".join(splitname[1:-1])
            # read csv with expression data into dataframe and get gene fold changes
            mag_fold1 = pd.read_csv(degpath + '/' + file1 + ".csv")
            print(mag_fold1.head())
            gene_ids = df['gene_id'].tolist()
            foldchanges = mag_fold1[mag_fold1['gene_id'].isin(gene_ids)]['log2FoldChange'].astype(float)
            # mag_fold1 = mag_fold1[['log2FoldChange']].astype(float)
            # print(mag_fold1.head())
            # add stuff into df
            # df['log2FoldChange'] = mag_fold1['log2FoldChange'].tolist()
            df['log2FoldChange'] = foldchanges
            df['dataset'] = file1
            df['overlap'] = file1 + ' (n:' + n + ')'
            # append dataframe to overall list for later visualization
            li.append(df)
        else:   # for intersection of genes between two experimental conditions, ie. overlapping group
            # read bed file into dataframe
            intersection = True
            splitname = filename.split('/')[-1].replace('.', '_').split("_")
            print(splitname)
            mid = int(len(splitname)/2)
            # get base names of each files, shorter if embryo which is unsexed
            if 'e' in temp:
                fullfile = "_".join(splitname[1:-1])
                file1 = "_".join(splitname[1:mid])
                file2 = "_".join(splitname[mid:-1])
            else:
                fullfile = "_".join(splitname[1:-1])
                file1 = "_".join(splitname[1:mid])
                file2 = "_".join(splitname[mid:-1])
            # read csv with expression data for first experimental condition into dataframe and get gene fold changes
            df = pd.read_csv(filename, sep="\t", names=['chr', 'start', 'end', 'gene_id', 'score', 'gene_name'])
            n = str(len(df.index))
            mag_fold1 = pd.read_csv(degpath + '/' + file1 + ".csv")
            gene_ids = df['gene_id'].tolist()
            foldchanges = mag_fold1[mag_fold1['gene_id'].isin(gene_ids)]['log2FoldChange'].astype(float)
            # add stuff to df
            df['log2FoldChange'] = foldchanges
            df['dataset'] = file1
            df['overlap'] = fullfile + ' (n:' + n + ')'
            # append dataframe for first set to overall list for later visualization
            li.append(df)
            # read csv with expression data for first experimental condition into dataframe and get gene fold changes
            df = pd.read_csv(filename, sep="\t", names=['chr', 'start', 'end', 'gene_id', 'score', 'gene_name'])
            n = str(len(df.index))
            mag_fold2 = pd.read_csv(degpath + '/' + file2 + ".csv")
            foldchanges = mag_fold2[mag_fold2['gene_id'].isin(gene_ids)]['log2FoldChange'].astype(float)
            # add stuff to df
            df['log2FoldChange'] = foldchanges
            df['dataset'] = file2
            df['overlap'] = fullfile + ' (n:' + n + ')'
            # append dataframe for first set to overall list for later visualization
            li.append(df)

    # concatenate list into single larger dataframe
    frame = pd.concat(li, axis=0, ignore_index=True)

    # generate empty plot 
    plt.figure(figsize=(9,5))
    
    if not intersection:
        unique_datasets = frame['dataset'].unique()

        if len(unique_datasets) >= 2:
            # determine up and down regulated gene sets
            done_up = False
            done_down = False
            upreg_set = unique_datasets[0]
            downreg_set = unique_datasets[1]
            for d in unique_datasets:
                if 'upreg' in d and not done_up:
                    upreg_set = d
                    done_up = True
                elif 'downreg' in d and not done_down:
                    downreg_set = d
                    done_down = False
                elif 'upreg' in d and done_up:
                    downreg_set = d
                    done_down = False
                elif 'downreg' in d and done_down:
                    upreg_set = d
                    done_up = True

            # set colors for different groups
            my_pal = {upreg_set: "#af8dc3", downreg_set: "#7fbf7b"}

    elif intersection:

        unique_datasets = frame['overlap'].unique()

        if len(unique_datasets) >= 2:
            # determine up and down regulated gene sets
            done_up = False
            done_down = False
            upreg_set = unique_datasets[0]
            downreg_set = unique_datasets[1]
            for d in unique_datasets:
                if 'upreg' in d and not done_up:
                    upreg_set = d
                    done_up = True
                elif 'downreg' in d and not done_down:
                    downreg_set = d
                    done_down = False
                elif 'upreg' in d and done_up:
                    downreg_set = d
                    done_down = False
                elif 'downreg' in d and done_down:
                    upreg_set = d
                    done_up = True
            # my_pal = {"clampi_e_upreg (n:456)": "#af8dc3", "clampi_e_downreg (n:377)":"#7fbf7b"}
            # my_pal = {"msl2i_A_m_upreg": "#af8dc3", "msl2i_A_m_downreg": "#7fbf7b"}

            # set colors for different groups
            my_pal = {upreg_set: "#af8dc3", downreg_set: "#7fbf7b"}

    # get all unique datasets (experimental conditions) from dataframe, frame
    unique_datasets = frame['dataset'].unique()
    # create combinations of size 2
    subsets3 = list(itertools.combinations(unique_datasets, 2))
    print(subsets3)
    # box_pairs_lst = []
    # for pair in subsets3:
    #     tup1,tup2 = pair
    #     ds1,chr1 = tup1
    #     ds2,chr2 = tup2
    #     if ds1 == ds2 or chr1 == chr2:
    #         box_pairs_lst.append(pair)
    box_pairs_lst = subsets3
    print('box pairs ', box_pairs_lst)

    # get all possible subsets of these combinations
    all_ss = []
    for L in range(1, len(box_pairs_lst)+1):
        for ss in itertools.combinations(box_pairs_lst, L):
            all_ss.append(list(ss))

    # sort subsets for largest to smallest by length
    all_ss_sorted = sorted(all_ss, key=len, reverse=True)

    print('all ss ', all_ss_sorted)

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
                plt.tight_layout(pad=2)
                plt.tick_params(axis='x', pad=17)
                # add sample size (n: _) under each of the violin groups
                nobs =  frame.groupby(['dataset']).apply(lambda x: 'n: {}'.format(len(x)))
                ymin,ymax = g.get_ylim()
                for ax in plt.gcf().axes:
                    for tick, label in enumerate(ax.get_xticklabels()):
                        ax_dat = label.get_text()
                        x_offset = 0
                        num = nobs[ax_dat]
                        ax.text(tick + x_offset, ymin - abs(0.0686*(ymax-ymin)), num,  #ymin - 0.25
                                horizontalalignment='center', size='medium', color='green', weight='semibold')
                # get title and save figure
                fig = g.get_figure()
                # if 'e' in splitname:
                #     figtit = "_".join(splitname[2:5])
                # else:
                #     figtit = "_".join(splitname[2:6])
                figtit = "__".join(frame['dataset'].unique())
                fig.savefig(savepath + '/overall_violin_'+figtit+'.png')
            elif intersection:
                # make violin plot of dataset (experimental condition) vs log2FoldChange
                print(my_pal)
                print(frame['overlap'].unique())
                plt.figure(figsize=(9,5))
                g = seaborn.violinplot(x="dataset", y="log2FoldChange", data=frame, hue="overlap", palette=my_pal)
                plt.legend(title='Groups with Same Genes', bbox_to_anchor=(1.05,1), loc=2)
                g.set_xticklabels(g.get_xticklabels(), rotation=15)
                # add significance bars
                add_stat_annotation(g, data=frame, x="dataset", y="log2FoldChange",
                                    box_pairs=blabla,
                                    test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
                plt.tight_layout(pad=2)
                # get title and save figure
                fig = g.get_figure()
                # if 'e' in splitname:
                #     figtit = "_".join(splitname[1:-1])
                # else:
                #     figtit = "_".join(splitname[1:-1])
                titlst = []
                for el in frame['overlap'].unique():
                    overlap_name = el.split(' ')[0]
                    titlst.append(overlap_name)
                figtit = "__".join(titlst)
                fig.savefig(savepath + '/overall_violin_'+figtit+'.png')
            break
        except ValueError:
            continue

    # plt.tight_layout(pad=2)

    # plt.tick_params(axis='x', pad=17)

    # # add sample size (n: _) under each of the violin groups
    # medians = frame.groupby(['dataset'])['log2FoldChange'].median() 
    # mins = frame.groupby(['dataset'])['log2FoldChange'].min() 
    # nobs =  frame.groupby(['dataset']).apply(lambda x: 'n: {}'.format(len(x)))
    # ymin,ymax = g.get_ylim()
    # for ax in plt.gcf().axes:
    #     for tick, label in enumerate(ax.get_xticklabels()):
    #         ax_dat = label.get_text()
    #         x_offset = 0
    #         min_val = mins[ax_dat]
    #         num = nobs[ax_dat]
    #         ax.text(tick + x_offset, ymin - abs(0.0686*(ymax-ymin)), num,  #ymin - 0.25
    #                 horizontalalignment='center', size='medium', color='green', weight='semibold')

    # get title and save figure
    # fig = g.get_figure()
    # if 'e' in splitname:
    #     figtit = "_".join(splitname[2:5])
    # else:
    #     figtit = "_".join(splitname[2:6])
    # fig.savefig('overall_violin_'+figtit+'.png')



def main():
    # create the command line parser
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--filepath", default='hello.bed', help="path to bed file to get basic gene data")
    parser.add_argument("-d", "--degpath", default='hello.csv', help="path to csv file to get deg data")
    parser.add_argument("-s", "--savepath", default='.', help="path to save figures")

    args = parser.parse_args()
    data = args.filepath
    degpath = args.degpath
    savepath = args.savepath
    
    # get list of bed files with genes sets to compare
    all_files = glob.glob(data + "/*.bed")
    print(savepath)

    if not os.path.isdir(savepath):
        os.makedirs(savepath)
        print('dir made')
	
    for idx,f in enumerate(all_files):
        print(f)
	#makeViolinPlot([f], degpath, savepath)
        remaining_files = all_files[idx+1:]
        for f2 in remaining_files:
            file_lst = [f,f2]
            print('file list: ',file_lst)
            makeViolinPlot(file_lst, degpath, savepath)
            # break
        # break


if __name__ == "__main__":
    main()
