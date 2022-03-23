# global_XA.py
# Purpose: plot global X v.s. autosomes
# Last mod. 02/22/2022

import pandas as pd
import numpy as np
import argparse
import glob
import seaborn
from statannot import add_stat_annotation
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import itertools
from scipy.stats import ttest_ind, mannwhitneyu, ks_2samp, anderson_ksamp
from statsmodels.stats.multitest import multipletests
from os import listdir
from os.path import isfile, join
import warnings
import os
warnings.filterwarnings('ignore')

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def makeGlobalPlot(filename, reg, savepath, gtf_path):
    li = []
    pth = filename.split('/')[:-1]
    pth = '/'.join(pth)
    df = pd.read_csv(filename)
    # expression_replicates = df[[caseid+'.1', caseid+'.2', caseid+'.3', caseid+'.4']]
    # avg_exp = np.sum(expression_replicates, axis=1) / 4
    # row_drop_genes_replicates = np.where(avg_exp <= 3)[0]
    # df2 = df.drop(df.index[row_drop_genes_replicates])
    if reg == 'down':
        row_drop_genes_foldchange = np.where(df['log2FoldChange'] > 0 )[0]
        df2 = df.drop(df.index[row_drop_genes_foldchange])
    elif reg == 'up':
        row_drop_genes_foldchange = np.where(df['log2FoldChange'] < 0 )[0]
        df2 = df.drop(df.index[row_drop_genes_foldchange])
    else:
        df2 = df
    # elif reg == 'all':
    #     pass        
    if 'chrom' not in df2.columns:
        # print(df2[0:2])
        gene_ids = df2['gene_id'].tolist()
        genes = pd.read_csv(gtf_path, sep="\t", names=['chr', 'start', 'end', 'gene_id', 'score', 'gene_name'])
        # print(genes[0:10])
        # print(gene_ids[0:10])
        # print(genes['gene_id'].isin(gene_ids)[0:10])
        # print(genes[genes['gene_id'].isin(gene_ids)][0:10])
        # overlaps = genes[genes['gene_id'].isin(gene_ids)]
        # chromosomes = overlaps['chr'].tolist()
        chromosomes = []
        for gt in gene_ids:
            # print(als)
            rowchr = genes.loc[genes['gene_id'] == gt]
            # print(rowchr)
            actchr = rowchr['chr'].values[0]
            # print(actchr)
            chromosomes.append(actchr)
        foldchanges = df2['log2FoldChange'].tolist()
        pvals = df2['padj'].tolist()
        new_df3 = pd.DataFrame({'chrom':chromosomes,'log2FoldChange':foldchanges,'pvalue':pvals})
        new_df3 = new_df3[new_df3['chrom'].isin(['2L','2R','3L','3R','4','X'])]
        li.append(new_df3)
    else:
        chromosomes = df2['chrom'].tolist()
        foldchanges = df2['log2FoldChange'].tolist()
        pvals = df2['padj'].tolist()
        new_df3 = pd.DataFrame({'chrom':chromosomes,'log2FoldChange':foldchanges,'pvalue':pvals})
        new_df3 = new_df3[new_df3['chrom'].isin(['2L','2R','3L','3R','4','X'])]
        li.append(new_df3)
    df4_A = new_df3.copy()
    df4_A['chrom'] = df4_A['chrom'].map({'2L':'A', '2R':'A', '3L':'A', '3R':'A', '4':'A', 'X':'X'})
    row_drop_genes_X = np.where(new_df3['chrom'] == 'X')[0]
    df4_A = df4_A.drop(df4_A.index[row_drop_genes_X])
    li.append(df4_A)

    frame = pd.concat(li, axis=0, ignore_index=True)

    frame = frame.sort_values('chrom')
    # print(frame.tail(15))
    print('about to make fig')
    # plt.figure(figsize=(8,5))
    f = plt.figure(figsize=(8,5))
    ax = f.add_subplot(111)
    seaborn.boxplot(x="chrom", y="log2FoldChange", data=frame, notch=True)
    print('made box plot')
    #cpalette = seaborn.color_palette("Reds_r")[:-1] #or Reds_r palette YlOrBr_r #other color


    # seaborn.swarmplot(x="chrom", y="log2FoldChange", hue="pvalue", data=frame, marker='d', palette="Reds_r", alpha=0.9, size=2) #this line is taking too long
    

    # g2 = seaborn.swarmplot(x="chrom", y="log2FoldChange", hue="pvalue", data=frame, marker='d', palette="Reds_r", alpha=0.7, size=2)
    # g2 = seaborn.swarmplot(x="chrom", y="log2FoldChange", hue="pvalue", marker='d', data=frame, palette=cpalette, alpha=0.6)
    # g2 = seaborn.swarmplot(x="chrom", y="log2FoldChange", hue="pvalue", marker='d', data=frame, color='red', alpha=0.6)
    print('added swarm plot on top')
    # g2.legend_.remove()
    # ax.get_legend().set_visible(False)
    ax.set(xlabel='Chromosomes', ylabel='log2FoldChange')
    plt.tick_params(axis='x', pad=17)

    plt.tight_layout()

    condition = filename.split('/')[-1].split('.')[0]
    plt.title('Global X vs A Plot for ' + condition + ' ' + reg)

    nobs =  frame.groupby(['chrom']).apply(lambda x: 'n: {}'.format(len(x)))
    ymin,ymax = ax.get_ylim()
    for axc in plt.gcf().axes:
        for tick, label in enumerate(axc.get_xticklabels()):
            ax_dat = label.get_text()
            x_offset = 0
            num = nobs[ax_dat]
            axc.text(tick + x_offset, ymin - abs(0.0486*(ymax-ymin)), num,
                    horizontalalignment='center', size='small', color='green', weight='semibold')

    #plt.show()
    # fig = g.get_figure()
    # pth = filename.split('/')[:-1]
    # pth = '/'.join(pth)
    fn = filename.split('/')[-1].split(".")[0]
    plt.savefig(savepath+'/globalXA_'+fn+'_'+reg+'.png')

    plt.clf()

    chromosomes = list(frame['chrom'].unique())
    if 'X' in chromosomes:
        chromosomes.insert(0, chromosomes.pop(chromosomes.index('X')))
    subsets = list(itertools.combinations(chromosomes, 2))

    comparisons_lst = []
    ttest_lst = []
    mannwhitney_lst = []
    kstest_lst = []
    anderson_teststat_lst = []
    anderson_critval_lst = []

    for chr1,chr2 in subsets:

        cat1 = frame[frame['chrom']==chr1]
        cat2 = frame[frame['chrom']==chr2]
        ttest_pval = ttest_ind(cat1['log2FoldChange'], cat2['log2FoldChange'])[1]
        mannwhitney_pval = mannwhitneyu(cat1['log2FoldChange'], cat2['log2FoldChange'])[1]
        kstest_pval = ks_2samp(cat1['log2FoldChange'], cat2['log2FoldChange'])[1]
        anderson_test = anderson_ksamp([cat1['log2FoldChange'], cat2['log2FoldChange']])
        anderson_teststat = anderson_test[0]
        anderson_critval = anderson_test[1][2]
        comp = chr1 + ' and ' + chr2
        if 'X' in comp and '4' in comp or 'X' in comp and 'A' in comp:
            comp = comp + '  '
        elif 'X' in comp or '4' in comp or 'A' in comp:
            comp = comp + ' '
        comparisons_lst.append(comp)
        # ttest_lst.append(float(str(ttest_pval)[:6]))
        # mannwhitney_lst.append(float(str(mannwhitney_pval)[:6]))
        # kstest_lst.append(float(str(kstest_pval)[:6]))
        ttest_lst.append(ttest_pval)
        mannwhitney_lst.append(mannwhitney_pval)
        kstest_lst.append(kstest_pval) 
        anderson_teststat_lst.append(anderson_teststat)
        anderson_critval_lst.append(anderson_critval)

    print(anderson_teststat_lst)
    print(anderson_critval_lst)
    #output critical value and test statistic, make sure test statistic is greater than critical value, if test stat is greater than crit value then reject null
    anderson_lst = []
    anderson_lst_binary = []
    for i in range(len(anderson_teststat_lst)):
        teststat = anderson_teststat_lst[i]
        critval = anderson_critval_lst[i]
        if teststat > critval:
            anderson_lst.append('reject null')
            anderson_lst_binary.append(1)
        else:
            anderson_lst.append('accept null')
            anderson_lst_binary.append(0)
    print(anderson_lst)
    print(ttest_lst)
    # bonferroni multiply by number of tests to correct for alpha cutoff
    ttest_lst = multipletests(ttest_lst, method='fdr_bh')[1].tolist() #bonferroni, hs, holm, fdr_bh # ***people usually use benjamini hoshberg bh***
    mannwhitney_lst = multipletests(mannwhitney_lst, method='fdr_bh')[1].tolist()
    kstest_lst = multipletests(kstest_lst, method='fdr_bh')[1].tolist()
    print(ttest_lst)

    data = {'comparisons': comparisons_lst,
            'T-test p-values    ': ttest_lst,
            'Mann Whitney p-vals': mannwhitney_lst,
            'KS-test p-values   ': kstest_lst,
            'AD Test Statistic  ': anderson_teststat_lst,
            'AD Critical Value  ': anderson_critval_lst,
            'AD Significance?   ': anderson_lst}
    stats_df = pd.DataFrame(data)
    print(stats_df)
    stats_df.to_csv(savepath+'/globalXA_'+fn+'_stats_'+reg+'.txt', sep='\t')
    

    stats_df2 = pd.DataFrame({'t-test p-vals': ttest_lst,
                            'mann-whitney p-vals': mannwhitney_lst,
                            'ks-test p-vals': kstest_lst}, 
                            index = comparisons_lst)
    stats_df2.rename_axis('comparisons')
    f = plt.figure(figsize=(7,5))
    ax = f.add_subplot(111)
    # cmap reds -> darkest red is most significant (closest to 0)
    pval_heatmap = seaborn.heatmap(stats_df2, annot=True, cmap = "Reds_r", vmin=0, vmax=1) #or rocket_r or coolwarm_r
    fig2 = pval_heatmap.get_figure()
    fig2.tight_layout()
    #plt.show(fig2)
    fig2.savefig(savepath+'/globalXA_'+fn+'_stats_'+reg+'.png')


    # stats_df3 = pd.DataFrame({'AD Significance?': anderson_lst_binary}, 
    #                     index = comparisons_lst)
    # # fig, ax = plt.subplots(1,2)
    # fig, (ax,ax2) = plt.subplots(ncols=2)
    # fig.subplots_adjust(wspace=0.01)
    # seaborn.heatmap(stats_df2, ax=ax,annot=True, cmap = "Reds_r")
    # seaborn.heatmap(stats_df3, ax=ax2, cmap="Greens", yticklabels=False, cbar=False, cbar_ax=None, square=True)
    # plt.show()


def main():
    # create the command line parser
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--filepath", default='C:/Users/isaac/OneDrive/Documents/BIOL_1950_Larschan_Lab', help="path to csv file to get RNA expression data")
    parser.add_argument("-s", "--savepath", default='.', help="path to save figures")
    parser.add_argument("-g", "--gtfpath", default='.', help="path to gtf bed file")

    args = parser.parse_args()
    filepath = args.filepath
    savepath = args.savepath
    gtf_path = args.gtfpath

    onlyfiles = [f for f in listdir(filepath) if isfile(join(filepath, f))]
    files_lst = []

    if not os.path.isdir(savepath):
        os.mkdir(savepath)

    combined = False
    for fil in onlyfiles:
        if 'combined' in fil:
            combined = True
        if fil.lower().endswith('.csv'):
            print(fil)
            files_lst.append(fil)

    for filename in files_lst:
        print(filename)
        if combined:
            regulation = ['all', 'up', 'down']
            for r in regulation:
                print(r)
                makeGlobalPlot(filepath+'/'+filename, r, savepath, gtf_path)
        else:
            if 'upreg' in filename:
                makeGlobalPlot(filepath+'/'+filename, 'up', savepath, gtf_path)
            elif 'downreg' in filename:
                makeGlobalPlot(filepath+'/'+filename, 'down', savepath, gtf_path)
            else:
                makeGlobalPlot(filepath+'/'+filename, 'all', savepath, gtf_path)

if __name__ == "__main__":
    main()