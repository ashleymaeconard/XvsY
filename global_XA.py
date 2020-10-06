import pandas as pd
import numpy as np
import argparse
import glob
import seaborn
from statannot import add_stat_annotation
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import itertools
from scipy.stats import ttest_ind, mannwhitneyu, ks_2samp

import warnings
warnings.filterwarnings('ignore')

# truncate colormap range to get rid of colors that are too light or too dark
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def main(filename, caseid, reg):
    li = []
    # get file path
    pth = filename.split('/')[:-1]
    pth = '/'.join(pth)
    # read in expression data from DESeq2 from CSV
    df = pd.read_csv(filename)
    # filter through replicates to remove those with mean counts less than 3
    expression_replicates = df[[caseid+'.1', caseid+'.2', caseid+'.3', caseid+'.4']]
    avg_exp = np.sum(expression_replicates, axis=1) / 4
    row_drop_genes_replicates = np.where(avg_exp <= 3)[0]
    df2 = df.drop(df.index[row_drop_genes_replicates])
    # filter for downregulated genes, upregulated genes, or do not filter (keep all)
    if reg == 'down':
        row_drop_genes_foldchange = np.where(df2['log2FoldChange'] > 0 )[0]
        df2 = df2.drop(df2.index[row_drop_genes_foldchange])
    elif reg == 'up':
        row_drop_genes_foldchange = np.where(df2['log2FoldChange'] < 0 )[0]
        df2 = df2.drop(df2.index[row_drop_genes_foldchange])
    elif reg == 'all':
        pass
    # if chromosomes not found in dataframe, search through gtf_genes.bed file to get this information for the
    # corresponding FlyBase IDs
    if 'chrom' not in df2.columns:
        gene_ids = df2['gene_id'].tolist()
        genes = pd.read_csv(pth+"/gtf_genes.bed", sep="\t", names=['chr', 'start', 'end', 'gene_id', 'score', 'gene_name'])
        overlaps = genes[genes['gene_id'].isin(gene_ids)]
        chromosomes = overlaps['chr'].tolist()
        foldchanges = df2['log2FoldChange'].tolist()
        pvals = df2['padj'].tolist()
        # construct new dataframe with chrom, fold change, and pval
        new_df3 = pd.DataFrame({'chrom':chromosomes,'log2FoldChange':foldchanges,'pvalue':pvals})
        # filter in case other types of DNA (weird artifact)
        new_df3 = new_df3[new_df3['chrom'].isin(['2L','2R','3L','3R','4','X'])]
        # append for later visualization
        li.append(new_df3)
    else: # if dataframe has chromosome info
        chromosomes = df2['chrom'].tolist()
        foldchanges = df2['log2FoldChange'].tolist()
        pvals = df2['padj'].tolist()
        # construct new dataframe with chrom, fold change, and pval
        new_df3 = pd.DataFrame({'chrom':chromosomes,'log2FoldChange':foldchanges,'pvalue':pvals})
        # filter in case other types of DNA (weird artifact)
        new_df3 = new_df3[new_df3['chrom'].isin(['2L','2R','3L','3R','4','X'])]
        # append for later visualization
        li.append(new_df3)
    # create copy of previous dataframe
    df4_A = new_df3.copy()
    # change all chromosome labels to either X or A (autosome) and then drop the X
    df4_A['chrom'] = df4_A['chrom'].map({'2L':'A', '2R':'A', '3L':'A', '3R':'A', '4':'A', 'X':'X'})
    row_drop_genes_X = np.where(new_df3['chrom'] == 'X' )[0]
    df4_A = df4_A.drop(df4_A.index[row_drop_genes_X])
    # append these autosomes to the list
    li.append(df4_A)

    # concatenate list into single larger dataframe
    frame = pd.concat(li, axis=0, ignore_index=True)

    # sort by chromosome to have A and X at the end of plot and 2L at beginning
    frame = frame.sort_values('chrom')

    # generate seaborn notched box plot for chrom vs fold change
    plt.figure(figsize=(8,5))
    g = seaborn.boxplot(x="chrom", y="log2FoldChange", data=frame, notch=True)

    # overlay swarmplot with pvalue of individual genes
    cpalette = seaborn.color_palette("YlOrBr_r")[:-1] #or Reds_r palette YlOrBr_r
    g2 = seaborn.swarmplot(x="chrom", y="log2FoldChange", hue="pvalue", data=frame, palette=cpalette, alpha=0.6)
    g2.legend_.remove()
    plt.tick_params(axis='x', pad=17)

    plt.tight_layout()

    # add sample size (n: _) under each of the violin groups
    nobs =  frame.groupby(['chrom']).apply(lambda x: 'n: {}'.format(len(x)))
    ymin,ymax = g.get_ylim()
    for ax in plt.gcf().axes:
        for tick, label in enumerate(ax.get_xticklabels()):
            ax_dat = label.get_text()
            x_offset = 0
            num = nobs[ax_dat]
            ax.text(tick + x_offset, ymin - abs(0.0486*(ymax-ymin)), num,
                    horizontalalignment='center', size='small', color='green', weight='semibold')

    plt.show()
    fig = g.get_figure()
    # pth = filename.split('/')[:-1]
    # pth = '/'.join(pth)

    # get filename and save figure
    fn = filename.split('/')[-1].split(".")[0]
    # fig.savefig(pth+'/globalXA_'+fn+'_'+reg+'.png')

    # get list of chromosomes and move X to the front
    chromosomes = list(frame['chrom'].unique())
    if 'X' in chromosomes:
        chromosomes.insert(0, chromosomes.pop(chromosomes.index('X')))
    # create all pairs of combinations of chromsomes
    subsets = list(itertools.combinations(chromosomes, 2))

    # initialize lists for statistical tests
    comparisons_lst = []
    ttest_lst = []
    mannwhitney_lst = []
    kstest_lst = []

    # run statistical tests on each combination of chromosomes
    for chr1,chr2 in subsets:

        cat1 = frame[frame['chrom']==chr1]
        cat2 = frame[frame['chrom']==chr2]
        ttest_pval = ttest_ind(cat1['log2FoldChange'], cat2['log2FoldChange'])[1]
        mannwhitney_pval = mannwhitneyu(cat1['log2FoldChange'], cat2['log2FoldChange'])[1]
        kstest_pval = ks_2samp(cat1['log2FoldChange'], cat2['log2FoldChange'])[1]
        comp = chr1 + ' vs ' + chr2
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

    # create dataframe with these pvalue results and save as txt file
    data = {'comparisons': comparisons_lst,
            'T-test p-values    ': ttest_lst,
            'Mann Whitney p-vals': mannwhitney_lst,
            'KS-test p-values   ': kstest_lst}
    stats_df = pd.DataFrame(data)
    # stats_df.to_csv(pth+'/globalXA_'+fn+'_stats_'+reg+'.txt', sep='\t')

    # make second dataframe to create heatmap of pvalues
    stats_df2 = pd.DataFrame({'t-test p-vals': ttest_lst,
                            'mann-whitney p-vals': mannwhitney_lst,
                            'ks-test p-vals': kstest_lst}, 
                            index = comparisons_lst)
    stats_df2.rename_axis('comparisons')

    pval_heatmap = seaborn.heatmap(stats_df2, annot=True, cmap = "coolwarm_r") #or rocket_r
    fig2 = pval_heatmap.get_figure()
    plt.show(fig2)
    # fig2.savefig(pth+'/globalXA_'+fn+'_stats_'+reg+'.png')


if __name__ == "__main__":
    # create the command line parser
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--file", default='C:/Users/isaac/OneDrive/Documents/BIOL_1950_Larschan_Lab', help="path to csv file to get RNA expression data")
    parser.add_argument("-c", "--caseid", default='C4e', help="experiment data case ID (ex C4e)")

    args = parser.parse_args()
    filename = args.file
    caseid = args.caseid

    # run for down, up and all regulated genes
    regulation = ['down', 'up', 'all']
    for r in regulation:
        main(filename, caseid, r)