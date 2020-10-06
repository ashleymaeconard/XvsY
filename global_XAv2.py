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
from os import listdir
from os.path import isfile, join
import warnings
warnings.filterwarnings('ignore')

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def main(filename, caseid, reg):
    li = []
    pth = filename.split('/')[:-1]
    pth = '/'.join(pth)
    df = pd.read_csv(filename)
    expression_replicates = df[[caseid+'.1', caseid+'.2', caseid+'.3', caseid+'.4']]
    avg_exp = np.sum(expression_replicates, axis=1) / 4
    row_drop_genes_replicates = np.where(avg_exp <= 3)[0]
    df2 = df.drop(df.index[row_drop_genes_replicates])
    if reg == 'down':
        row_drop_genes_foldchange = np.where(df2['log2FoldChange'] > 0 )[0]
        df2 = df2.drop(df2.index[row_drop_genes_foldchange])
    elif reg == 'up':
        row_drop_genes_foldchange = np.where(df2['log2FoldChange'] < 0 )[0]
        df2 = df2.drop(df2.index[row_drop_genes_foldchange])
    elif reg == 'all':
        pass        
    if 'chrom' not in df2.columns:
        gene_ids = df2['gene_id'].tolist()
        genes = pd.read_csv(pth+"/gtf_genes.bed", sep="\t", names=['chr', 'start', 'end', 'gene_id', 'score', 'gene_name'])
        overlaps = genes[genes['gene_id'].isin(gene_ids)]
        chromosomes = overlaps['chr'].tolist()
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
    row_drop_genes_X = np.where(new_df3['chrom'] == 'X' )[0]
    df4_A = df4_A.drop(df4_A.index[row_drop_genes_X])
    li.append(df4_A)

    frame = pd.concat(li, axis=0, ignore_index=True)

    frame = frame.sort_values('chrom')

    f = plt.figure(figsize=(8,5))
    ax = f.add_subplot(111)
    g = seaborn.boxplot(x="chrom", y="log2FoldChange", data=frame, notch=True)

    #cpalette = seaborn.color_palette("Reds_r")[:-1] #or Reds_r palette YlOrBr_r #other color
    g2 = seaborn.swarmplot(x="chrom", y="log2FoldChange", hue="pvalue", data=frame, marker='d', palette="Reds_r", alpha=0.9, size=2)
    # g2 = seaborn.swarmplot(x="chrom", y="log2FoldChange", hue="pvalue", data=frame, marker='d', palette="Reds_r", alpha=0.7, size=2)
    # g2 = seaborn.swarmplot(x="chrom", y="log2FoldChange", hue="pvalue", marker='d', data=frame, palette=cpalette, alpha=0.6)
    # g2 = seaborn.swarmplot(x="chrom", y="log2FoldChange", hue="pvalue", marker='d', data=frame, color='red', alpha=0.6)
    g2.legend_.remove()
    #ax.get_legend().set_visible(False)
    #ax.set(xlabel='Chromosomes', ylabel='log2FoldChange')
    plt.tick_params(axis='x', pad=17)

    plt.tight_layout()

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
    fn = filename.split('/')[-1].split(".")[0]
    # fig.savefig(pth+'/globalXA_'+fn+'_'+reg+'.png')

    chromosomes = list(frame['chrom'].unique())
    if 'X' in chromosomes:
        chromosomes.insert(0, chromosomes.pop(chromosomes.index('X')))
    subsets = list(itertools.combinations(chromosomes, 2))

    comparisons_lst = []
    ttest_lst = []
    mannwhitney_lst = []
    kstest_lst = []

    for chr1,chr2 in subsets:

        cat1 = frame[frame['chrom']==chr1]
        cat2 = frame[frame['chrom']==chr2]
        ttest_pval = ttest_ind(cat1['log2FoldChange'], cat2['log2FoldChange'])[1]
        mannwhitney_pval = mannwhitneyu(cat1['log2FoldChange'], cat2['log2FoldChange'])[1]
        kstest_pval = ks_2samp(cat1['log2FoldChange'], cat2['log2FoldChange'])[1]
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

    data = {'comparisons': comparisons_lst,
            'T-test p-values    ': ttest_lst,
            'Mann Whitney p-vals': mannwhitney_lst,
            'KS-test p-values   ': kstest_lst}
    stats_df = pd.DataFrame(data)
    # stats_df.to_csv(pth+'/globalXA_'+fn+'_stats_'+reg+'.txt', sep='\t')


    stats_df2 = pd.DataFrame({'t-test p-vals': ttest_lst,
                            'mann-whitney p-vals': mannwhitney_lst,
                            'ks-test p-vals': kstest_lst}, 
                            index = comparisons_lst)
    stats_df2.rename_axis('comparisons')
    # cmap reds -> darkest red is most significant (closest to 0)
    pval_heatmap = seaborn.heatmap(stats_df2, annot=True, cmap = "Reds_r") #or rocket_r or coolwarm_r
    fig2 = pval_heatmap.get_figure()
    plt.show(fig2)
    # fig2.savefig(pth+'/globalXA_'+fn+'_stats_'+reg+'.png')


if __name__ == "__main__":
    # create the command line parser
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--filepath", default='C:/Users/isaac/OneDrive/Documents/BIOL_1950_Larschan_Lab', help="path to csv file to get RNA expression data")
    parser.add_argument("-m", "--metadata", default='C:/Users/isaac/OneDrive/Documents/BIOL_1950_Larschan_Lab', help="path to csv file to case ID and condition info")
    # parser.add_argument("-c", "--caseid", default='C4e', help="experiment data case ID (ex C4e)")

    args = parser.parse_args()
    filepath = args.filepath
    metafile = args.metadata
    # caseid = args.caseid

    metadata = pd.read_csv(metafile)
    all_cases = metadata['ID'].tolist()
    all_cond = metadata['condition'].tolist()

    all_cond_split = []
    for condition in all_cond:
        condition = condition.split("_")[1:]
        condition = "_".join(condition)
        all_cond_split.append(condition)

    caseid_lst = []

    onlyfiles = [f for f in listdir(filepath) if isfile(join(filepath, f))]
    files_lst = []

    combined = False
    for fil in onlyfiles:
        if 'combined' in fil:
            combined = True
        if fil.lower().endswith('.csv'):
            print(fil)
            split_fil = fil.split("_")[1:]
            split_fil = "_".join(split_fil)
            split_fil = split_fil.split(".")[0]
            ind = all_cond_split.index(split_fil)
            caseid_exp = all_cases[ind]
            caseid_lst.append(caseid_exp)
            files_lst.append(fil)

    for filename,caseid in zip(files_lst,caseid_lst):
        if combined:
            regulation = ['down', 'up', 'all']
            for r in regulation:
                main(filepath+'/'+filename, caseid, r)
        else:
            if 'upreg' in filename:
                main(filepath+'/'+filename, caseid, 'up')
            elif 'downreg' in filename:
                main(filepath+'/'+filename, caseid, 'down')