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

# Generate table with bold column headings and colored rows for FIMO results
def render_fimo_table(data, col_width=3.0, row_height=0.625, font_size=14,
                            header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                            bbox=[0, 0, 1, 1], header_columns=0,
                            ax=None, **kwargs):
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=size)
        # ax = plt.subplot(111)
        ax.axis('off')

    fimo_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)

    fimo_table.auto_set_font_size(False)
    fimo_table.set_fontsize(font_size)

    for k, cell in  six.iteritems(fimo_table._cells):
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        else:
            cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
    fig = ax.get_figure()
    fig.savefig(path+'\summary_table2.png')
    return ax


def main(path):

    # Gets all subdirectories in which MEME/FIMO results are found
    subdirs = [os.path.abspath(name) for name in os.listdir(path) if os.path.isdir(name)]

    # If there are no subdirectors (direct file input)
    if len(subdirs) == 0:

        # Get FIMO results (fimo.tsv)
        data = path + '/fimo.tsv'

        # Read TSV file into DataFrame and filter data
        df = pd.read_csv(data, sep="\t", names=['motif_id', 'motif_alt_id', 'sequence_name', 'start', 'stop', 'strand', 'score', 'p-value', 'q-value', 'matched_sequence'])
        df = df[['motif_id','sequence_name']]
        temp = df['sequence_name'].str.split("::", n = 1, expand=True)
        df['gene_name']=temp[0]
        df.drop(columns=['sequence_name'], inplace=True)
        df = df.iloc[:-3]
        df = df.iloc[1:]
        # Drop repeats in df where same gene has the same motif to get unique occurences
        df.drop_duplicates(subset=['motif_id','gene_name'], keep='first', inplace=True)
        df = df.sort_values(by = 'gene_name')
        df = df.reset_index()
        df.drop(columns=['index'], inplace=True)
        # Drop 'strange' motif that comes out of CLAMP PWM
        df = df[~df['motif_id'].str.contains("GTGTGT")]
        print(df)

        # Get list of unique motifs
        motifs = df['motif_id'].unique()
        print(motifs)

        # Get list of unique genes
        genes = df['gene_name'].unique()
        print(genes)
        num_genes = len(genes)

        lst_motifs=[]
        lst_percent=[]

        # Calculate percent occurence of each motif
        for m in motifs:
            df2 = df[df['motif_id'] == m]
            matches = df2.shape[0]
            percent = matches/num_genes
            lst_motifs.append(m)
            lst_percent.append(percent)

        # Zip results into summary dataframe
        summary_df = pd.DataFrame(list(zip(lst_motifs,lst_percent)), columns=['Motif','% Matches'])
        print(summary_df)

        # Generate basic figure of results
        ax1 = plt.subplot(111,frame_on=False)
        ax1.xaxis.set_visible(False)
        ax1.yaxis.set_visible(False)
        table(ax1,summary_df, loc='center')
        fig1 = ax1.get_figure()
        fig1.tight_layout()
        fig1.savefig(path+'\summary_table.png', bbox_inches='tight')

        # Generate pretty figure of results
        render_fimo_table(summary_df, header_columns=0, col_width=3.5)

    # If subdirectories exist
    elif len(subdirs) > 0:

        # For each gene cluster
        for cluster in subdirs:
        
            data = cluster + '/MEME/fimo.tsv'

            # Read TSV file into DataFrame and filter data
            df = pd.read_csv(data, sep="\t", names=['motif_id', 'motif_alt_id', 'sequence_name', 'start', 'stop', 'strand', 'score', 'p-value', 'q-value', 'matched_sequence'])
            df = df[['motif_id','sequence_name']]
            temp = df['sequence_name'].str.split("::", n = 1, expand=True)
            df['gene_name']=temp[0]
            df.drop(columns=['sequence_name'], inplace=True)
            df = df.iloc[:-3]
            df = df.iloc[1:]
            # Drop repeats in df where same gene has the same motif to get unique occurences
            df.drop_duplicates(subset=['motif_id','gene_name'], keep='first', inplace=True)
            df = df.sort_values(by = 'gene_name')
            df = df.reset_index()
            df.drop(columns=['index'], inplace=True)
            # Drop 'strange' motif that comes out of CLAMP PWM
            df = df[~df['motif_id'].str.contains("GTGTGT")]
            print(df)

            # Get list of unique motifs
            motifs = df['motif_id'].unique()
            print(motifs)

            # Get list of unique genes
            genes = df['gene_name'].unique()
            print(genes)
            num_genes = len(genes)

            lst_motifs=[]
            lst_percent=[]

            # Calculate percent occurence of each motif
            for m in motifs:
                df2 = df[df['motif_id'] == m]
                matches = df2.shape[0]
                percent = float(str((matches/num_genes) * 100)[:6])
                lst_motifs.append(m)
                lst_percent.append(percent)

            # Zip results into summary dataframe and sort by motif to have CLAMP first
            summary_df = pd.DataFrame(list(zip(lst_motifs,lst_percent)), columns=['Motif','% Matches'])
            summary_df = summary_df.sort_values(by = 'Motif', ascending=False)
            print(summary_df)

            # Generate basic figure of results
            ax1 = plt.subplot(111,frame_on=False)
            ax1.xaxis.set_visible(False)
            ax1.yaxis.set_visible(False)
            table(ax1,summary_df, loc='center')
            fig1 = ax1.get_figure()
            fig1.tight_layout()
            fig1.savefig(cluster+'/MEME/summary_table.png', bbox_inches='tight')

            # Generate pretty figure of results
            render_fimo_table(summary_df, header_columns=0, col_width=3.5)

if __name__ == "__main__":
    # create the command line parser
    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--path", default='.', help="path to meme results")

    args = parser.parse_args()
    path = args.path

    if path == '.':
        path = os.getcwd()

    # run main function to generate summary figures
    main(path)