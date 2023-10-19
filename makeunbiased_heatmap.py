'''
Filename: Create_unbiased_heatmap.py
Author: Joel Klein
Date: 2022-07--03
Description: This script creates a heatmap from a deseq2 output file and an annotation file with gene names and
creates a heatmap of the log2foldchange values.
License: CC-BY 4.0 International
Contact: joel.klein@wur.nl
'''

# Import packages
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


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

def create_heatmap(df):
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
    df.set_index('Names', inplace=True)
    # Set the color palette
    # Create the heatmap using Seaborn, set the z-scores
    #ax = sns.heatmap(df, cmap='RdBu', center=0, robust=True, cbar_kws={'label': 'z-scores', 'orientation': 'horizontal'})
    # create the heatmap without z-scores
    ax = sns.heatmap(df, cmap='RdBu_r', yticklabels=True, cbar_kws={'label': 'log2 fold change', 'orientation': 'horizontal'}, center=0)  # need to get log2foldchange around 0
    # Create the heatmap using Seaborn, set the z-scores and add the dentrodgram
    # set labels on the x and y axis
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=16)
    ax.set_ylabel('Gene ID', fontsize=22)
    figure = ax.get_figure()
    return figure


def annotate_dataframe(df, annotation_df):
    """
    Annotate the dataframe with the parasponia gene names
    :param df: deseq2 dataframe
    :return: deseq2 dataframe with parasponia gene names
    """
    # Select columns and remove 'Parand' and '.1' from Gene_id
    df = df[['Gene_id', 'nsp2-9'] + [col for col in df.columns if col not in ['Gene_id', 'nsp2-9']]]
    df['Gene_id'] = df['Gene_id'].str[7:-2]
    # Read in annotation data and select relevant columns
    annotation_df = annotation_df[annotation_df['# feature'].str.contains('mRNA')][['locus_tag', 'name', 'symbol']]
    annotation_df.rename(columns={'locus_tag': 'Gene_id'}, inplace=True)
    # Merge dataframes on Gene_id and create Annotation column
    df = pd.merge(df, annotation_df, on='Gene_id', how='left', suffixes=('', ''))
    df['Annotation'] = df['symbol'].fillna('') + ' | ' + df['name'].fillna('')
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


# read in the data
df_lfc = pd.read_csv('/path/to/input/LFC_table.tsv', sep='\t')
annotation_df = pd.read_csv('/path/to/input/LFC_table.csv',
    sep=',')
figure_name = "Output_heatmap"  # a name for the figure and output table

# Assuming df is your input DataFrame
df_lfc_annot = annotate_dataframe(df_lfc)
print(df_lfc_annot.head())

df_annot_merged = merge_geneid_annotation(df_lfc_annot)
# df_annot_merged.to_csv('/Users/joel/Documents/Students/nick/RNA-seq/NSP2_NP-_LFC_annotated.tsv', sep='\t', index=False, header=True)

# filter dataframe based on log2foldchange
df_lfdc_filtered = filter_lfc(df_annot_merged, -6, 6) # adjust min and max log2foldchange here

# create the heatmap
df_lfdc_filtered_sorted = df_lfdc_filtered.sort_values(by='NSP2ox6')
# save df to tsv
# df_lfdc_filtered_sorted.to_csv('/Users/joel/Documents/Students/nick/RNA-seq/{}_annotated_filtered.tsv'.format(figure_name), sep='\t', index=False)

# create the heatmap figures
df_lfc_fig = create_heatmap(df_lfdc_filtered_sorted)

# Configure matplotlib parameters for PDF saving
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = 'Arial'  # Set font to Arial

# Get the Figure object from the AxesSubplot object
# Adjust the margins of the figure
df_lfc_fig.subplots_adjust(left=0.4)
# Adjust the size of the figure
df_lfc_fig.set_size_inches(20, 80) # with, height
#plt.show()

# Save the Seaborn figure as a PDF file
# save fig with name figure_name
df_lfc_fig.savefig('/path/to/outputfile/pathwayheatmaps_{}_heatmapLFC.pdf'.format(figure_name), bbox_inches='tight')
print('Figure saved as PDF file')
df_lfc_fig.clf()  # remove the figure



