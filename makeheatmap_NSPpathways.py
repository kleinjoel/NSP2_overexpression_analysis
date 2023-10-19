'''
Filename: Create_unbiased_heatmap.py
Author: Joel Klein
Date: 2022-07--03
Description: This script creates a heatmap from a deseq2 output file and an annotation file with gene names and
creates a heatmap of the log2foldchange values.
License: CC-BY 4.0 International
Contact: joel.klein@wur.nl
'''

# import packages
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import zscore
import numpy as np


def filter_lfc(df_in, lfc_low, lfchigh):
    """
    Filter dataframe based on log2foldchange
    :param df_in: dataframe
    :param lfc_low: lower log2foldchange (Eg -2)
    :param lfchigh: higher log2foldchange (Eg 2)
    :return: filtered dataframe """
    print("filtering dataframe, this can take while...")
    df_filtered = df_in.loc[(df_in['NSP2ox1'] < lfc_low) | (df_in['NSP2ox3'] < lfc_low) | (df_in['NSP2ox6'] < lfc_low) | (df_in['NSP2ox1'] > lfchigh) | (df_in['NSP2ox3'] > lfchigh) | (df_in['NSP2ox6'] > lfchigh)]
    print("Done filtering dataframe, found {} DEGS with LFC cutoffs {},{}".format(len(df_filtered), lfc_low, lfchigh))
    return df_filtered

def annotate_dataframe(df, annotfile):
    """
    Annotate the dataframe with the parasponia gene names
    :param df: deseq2 dataframe
    :return: deseq2 dataframe with parasponia gene names
    """
    # Select columns and remove 'Parand' and '.1' from Gene_id
    df = df[['Gene_id', 'nsp2-9'] + [col for col in df.columns if col not in ['Gene_id', 'nsp2-9']]]
    df['Gene_id'] = df['Gene_id'].str[7:-2]
    # Read in annotation data and select relevant columns
    annotation_df = pd.read_csv(annotfile, sep=',')
    annotation_df = annotation_df[annotation_df['# feature'].str.contains('mRNA')][['locus_tag', 'name', 'symbol']]
    annotation_df.rename(columns={'locus_tag': 'Gene_id'}, inplace=True)
    # Merge dataframes on Gene_id and create Annotation column
    df = pd.merge(df, annotation_df, on='Gene_id', how='left', suffixes=('', ''))
    df['Annotation'] = df['symbol'].fillna('') + ' ' + df['name'].fillna('')
    df = df.drop(columns=['symbol', 'name'])
    # show all headers
    # pd.set_option('display.max_columns', None)
    # print(df.head())
    # Reorder columns
    new_order = ['Annotation', 'Gene_id', 'nsp2-9', 'NSP2ox1', 'NSP2ox3', 'NSP2ox6']
    df = df.reindex(columns=new_order)
    return df

def merge_geneid_annotation(df):
    """
    Merge Gene_id and Annotation columns
    :param df: dataframe with Gene_id and Annotation columns and log2foldchange columns
    :return df: dataframe with merged Gene_id and Annotation columns
    """
    df_copy = df.copy()
    df_copy.loc[:, 'Names'] = df_copy['Annotation'] + ' '+ df_copy['Gene_id']
    df_copy = df_copy.drop(columns=['Annotation', 'Gene_id'])
    # put column names as the first column
    cols = df_copy.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df_copy = df_copy[cols]
    # Remove leading "| " from values in the 'Annotation' column
    df_copy['Names'] = df_copy['Names'].str.lstrip('| ')
    print(df_copy.columns)
    return df_copy


def create_heatmapLFC(dfLFC, column_name):
    """
    Create a heatmap from a dataframe
    :param df: log2fold change dataframe
    :return: heatmap figure object
    """
    # Display the dataframe
    # pd.set_option('display.max_columns', None)
    # pd.set_option('display.max_rows', None)
    # print(df)
    # Set the gene_id column as the index
    print("Column name used for naming the samples:", column_name)
    dfLFC.set_index(column_name, inplace=True)
    # Set the color palette
    # Create the heatmap using Seaborn, set the z-scores
    #ax = sns.heatmap(df, cmap='RdBu', center=0, robust=True, cbar_kws={'label': 'z-scores', 'orientation': 'horizontal'})
    # create the heatmap without z-scores
    axlfc = sns.heatmap(dfLFC, cmap='RdBu_r', yticklabels=True, cbar_kws={'label': 'Log2FoldChange', 'orientation': 'horizontal'}, center=0)  # need to get log2foldchange around 0
    # Create the heatmap using Seaborn, set the z-scores and add the dentrodgram
    # set labels on the x and y axis
    # Set y-tick labels
    axlfc.set_yticklabels(axlfc.get_yticklabels(), fontsize=11)  # Set fontsize to 14 here
    # Set x-tick labels
    axlfc.set_xticklabels(axlfc.get_xticklabels(), fontsize=12, fontstyle='italic')
    # Set x-axis label
    axlfc.set_ylabel('Gene ID', fontsize=20)
    figurelfc = axlfc.get_figure()
    return figurelfc


def prepare_dataframe(cssp_filepath, df_lfc):
    """
    Merge pathway annotations with expression data and rearrange columns.
    :param cssp_filepath: Path to the pathway annotations file.
    :param df_lfc: DataFrame with log2 fold change information.
    :return: A processed DataFrame ready for heatmap plotting.
    """
    # Read in pathway annotations
    cssp_df = pd.read_csv(cssp_filepath, sep='\t')
    # Merge DataFrames and manipulate columns
    pathway_df = (
        pd.merge(cssp_df, df_lfc, left_on='geneid', right_on='Gene_id', how='inner')
            .assign(Gene_id=lambda x: x['annotation'] + ' | ' + x['geneid'])
            .drop(columns=['annotation', 'geneid'])
    )
    # Rearrange the columns to put 'Gene_id' first
    # Remove 'Parand_'
    pathway_df['Gene_id'] = pathway_df['Gene_id'].str.replace('Parand_', '')
    # Remove '.1' or '.2' at the end
    pathway_df['Gene_id'] = pathway_df['Gene_id'].str.replace(r'\.\d$', '')
    print(pathway_df.head())
    merged_pathway_df = pathway_df[['Gene_id'] + [col for col in pathway_df.columns if col != 'Gene_id']]
    # Further rearrange columns to move 'nsp2-9' after 'Gene_id'
    cols = ['Gene_id', 'nsp2-9'] + [col for col in merged_pathway_df.columns if col not in ['Gene_id', 'nsp2-9']]
    merged_pathway_df = merged_pathway_df[cols]
    # Sort by 'NSP2ox6'
    merged_pathway_dfsorted = merged_pathway_df.sort_values(by='NSP2ox6')
    # Print column names for debugging
    # print(merged_pathway_dfsorted.columns)
    return merged_pathway_dfsorted

################################################  MAIN OF THE SCRIPT ###################################################
# read in the data set figure names
df_lfc = pd.read_csv('/your/path/input/Deseq2_LFC.tsv', sep='\t')
pahtwayid_path = "/your/path/input/Pathwaytoplot.tsv"
figure_name= "egCarotenoid_NSP2_NP-"
########## ------------------------------##########--------------------------------##########-----------------##########

processed_df = prepare_dataframe(pahtwayid_path, df_lfc)

# make two copy to avoid error witht the dropped gene id column
processed_df_lfc = processed_df.copy()

## Create the LFC heatmap; Adjust the margins of the figure and save a high quality version for the LFC
df_lfc_fig = create_heatmapLFC(processed_df_lfc, "Gene_id")
df_lfc_fig.subplots_adjust(left=0.4)
# Adjust the size of the figure
df_lfc_fig.set_size_inches(15, 25) # with, height
df_lfc_fig.savefig('/your/path/toouput/{}_heatmaplfc.pdf'.format(figure_name), bbox_inches='tight')
print('Figure LFC saved as PDF file')
df_lfc_fig.clf()  # remove the figure

# merge two figures togheter ###################
figure_name= "P450_network_NSP2_NP-_NlowPhighmerged"
pahtwayid_path = "/your/path/input/Pathwaytoplot.tsv""

df_lfc_NPlow_df = pd.read_csv('/your/path/input/Deseq2_LFCLFC_condition1.tsv', sep='\t')
df_lfc_NlowPhigh = pd.read_csv('/your/path/input/Deseq2_LFCLFC_condition2.tsv', sep='\t')

processed_df_NPLOW = prepare_dataframe(pahtwayid_path, df_lfc_NPlow_df)
processed_df_NlowPhigh = prepare_dataframe(pahtwayid_path, df_lfc_NlowPhigh)

merged_df = processed_df_NPLOW.merge(processed_df_NlowPhigh, on='Gene_id', how='inner')

# Insert a NaN column for separation
merged_df['separator'] = np.nan
# Reorder the DataFrame columns
ordered_cols = ['Gene_id', 'nsp2-9_x', 'NSP2ox1_x', 'NSP2ox3_x', 'NSP2ox6_x', 'separator', 'nsp2-9_y', 'NSP2ox1_y', 'NSP2ox3_y', 'NSP2ox6_y']
merged_df = merged_df[ordered_cols]

# Configure matplotlib parameters for PDF saving
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = 'Arial'  # Set font to Arial

## Create the LFC heatmap; Adjust the margins of the figure and save a high quality version for the LFC
df_lfc_mergedfig = create_heatmapLFC(merged_df, "Gene_id")
df_lfc_mergedfig.subplots_adjust(left=0.4)
# Adjust the size of the figure
df_lfc_mergedfig.set_size_inches(22, 60)  # Width, Height in landscape for A4
df_lfc_mergedfig.savefig('//your/path/toouput//{}_heatmaplfc.pdf'.format(figure_name), bbox_inches='tight')
print('Figure LFC saved as PDF file')
df_lfc_mergedfig.clf()  # remove the figure
