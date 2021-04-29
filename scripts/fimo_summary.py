import pandas as pd
import numpy as np
import argparse
import glob
import seaborn
from statannot import add_stat_annotation
import matplotlib.pyplot as plt
import os
from pandas.plotting import table
import six

import warnings
warnings.filterwarnings('ignore')

# create the command line parser
parser = argparse.ArgumentParser()

parser.add_argument("-p", "--path", default='.', help="path to fimo results")
parser.add_argument("-s", "--save", default='.', help="where to save summary")

args = parser.parse_args()
path = args.path
save = args.save

if path == '.':
    path = os.getcwd()
if path[-1] != "/":
    path = path + "/"

subdirs = [path+name for name in os.listdir(path) if os.path.isdir(path+name)]

if len(subdirs) == 0:

    data = path + '/fimo.tsv'

    df = pd.read_csv(data, sep="\t", names=['motif_id', 'motif_alt_id', 'sequence_name', 'start', 'stop', 'strand', 'score', 'p-value', 'q-value', 'matched_sequence'])
    df = df[['motif_id','sequence_name']]
    temp = df['sequence_name'].str.split("::", n = 1, expand=True)
    df['gene_name']=temp[0]
    df.drop(columns=['sequence_name'], inplace=True)
    df = df.iloc[:-3]
    df = df.iloc[1:]
    df.drop_duplicates(subset=['motif_id','gene_name'], keep='first', inplace=True)
    df = df.sort_values(by = 'gene_name')
    df = df.reset_index()
    df.drop(columns=['index'], inplace=True)

    motifs = df['motif_id'].unique()
    print(motifs)

    genes = df['gene_name'].unique()
    num_genes = len(genes)

    lst_motifs=[]
    lst_percent=[]

    for m in motifs:
        df2 = df[df['motif_id'] == m]
        matches = df2.shape[0]
        percent = matches/num_genes
        lst_motifs.append(m)
        lst_percent.append(percent)

    summary_df = pd.DataFrame(list(zip(lst_motifs,lst_percent)), columns=['Motif','% Matches'])
    print(summary_df)

    # ax1 = plt.subplot(111,frame_on=False)
    # ax1.xaxis.set_visible(False)
    # ax1.yaxis.set_visible(False)
    # table(ax1,summary_df, loc='center')
    # fig1 = ax1.get_figure()
    # fig1.tight_layout()
    # fig1.savefig(save+'\summary_table.png', bbox_inches='tight')

    summary_df = pd.DataFrame(list(zip(lst_motifs,lst_percent)), columns=['Motif','% Matches in Genes'])
    def render_mpl_table(data, col_width=3.0, row_height=0.625, font_size=14,
                        header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                        bbox=[0, 0, 1, 1], header_columns=0,
                        ax=None, **kwargs):
        if ax is None:
            size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
            fig, ax = plt.subplots(figsize=size)
            # ax = plt.subplot(111)
            ax.axis('off')

        mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)

        mpl_table.auto_set_font_size(False)
        mpl_table.set_fontsize(font_size)

        for k, cell in  six.iteritems(mpl_table._cells):
            cell.set_edgecolor(edge_color)
            if k[0] == 0 or k[1] < header_columns:
                cell.set_text_props(weight='bold', color='w')
                cell.set_facecolor(header_color)
            else:
                cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
        fig = ax.get_figure()
        fig.savefig(save+'\summary_table.png')
        return ax

    render_mpl_table(summary_df, header_columns=0, col_width=3.5)

if len(subdirs) > 0:

    for cluster in subdirs:
        print(cluster)
        data = cluster + '/MEME/fimo.txt'
        df = pd.read_csv(data, sep="\t", names=['motif_id', 'sequence_name', 'start', 'stop', 'strand', 'score', 'p-value', 'q-value', 'matched_sequence'])
        print(df)
        df = df[1:]
        if df.empty:
            continue
        df = df[['motif_id','sequence_name']]
        temp = df['sequence_name'].str.split("::", n = 1, expand=True)
        df['gene_name']=temp[0]
        df.drop(columns=['sequence_name'], inplace=True)
        df.drop_duplicates(subset=['motif_id','gene_name'], keep='first', inplace=True)
        df = df.sort_values(by = 'gene_name')
        df = df.reset_index()
        df.drop(columns=['index'], inplace=True)

        motifs = df['motif_id'].unique()
        print(motifs)

        genes = df['gene_name'].unique()
        num_genes = len(genes)

        lst_motifs=[]
        lst_percent=[]

        for m in motifs:
            df2 = df[df['motif_id'] == m]
            matches = df2.shape[0]
            percent = (matches/num_genes) * 100
            lst_motifs.append(m)
            lst_percent.append(percent)

        summary_df = pd.DataFrame(list(zip(lst_motifs,lst_percent)), columns=['Motif','% Matches'])
        summary_df = summary_df.sort_values(by = 'Motif',ascending=False)
        print(summary_df)
        print(cluster)

        if motifs == [] or summary_df.empty:
            continue

        # ax1 = plt.subplot(111,frame_on=False)
        # ax1.xaxis.set_visible(False)
        # ax1.yaxis.set_visible(False)
        # table(ax1,summary_df, loc='center')
        # fig1 = ax1.get_figure()
        # fig1.tight_layout()
        # fig1.savefig(cluster+'/MEME/summary_table.png', bbox_inches='tight')

        summary_df = pd.DataFrame(list(zip(lst_motifs,lst_percent)), columns=['Motif','% Matches in Genes'])
        summary_df = summary_df.sort_values(by = 'Motif',ascending=False)
        def render_mpl_table(data, col_width=3.0, row_height=0.625, font_size=14,
                            header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                            bbox=[0, 0, 1, 1], header_columns=0,
                            ax=None, **kwargs):
            if ax is None:
                size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
                fig, ax = plt.subplots(figsize=size)
                # ax = plt.subplot(111)
                ax.axis('off')

            mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)

            mpl_table.auto_set_font_size(False)
            mpl_table.set_fontsize(font_size)

            for k, cell in  six.iteritems(mpl_table._cells):
                cell.set_edgecolor(edge_color)
                if k[0] == 0 or k[1] < header_columns:
                    cell.set_text_props(weight='bold', color='w')
                    cell.set_facecolor(header_color)
                else:
                    cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
            fig = ax.get_figure()
            fig.savefig(cluster+'/MEME/summary_table.png')
            return ax

        render_mpl_table(summary_df, header_columns=0, col_width=3.5)