# NSP2 overexpression analysis
Scripts for RNA-seq analysis and heatmap generation of NSP2 over-expression lines.

# Deseq2 analysis
Script name= Deseq2analysis.R 

This R script is designed to analyze Kallisto mapping output tables and perform a DESeq2 analysis on the two nutrient conditions separately. It uses two empty vector control lines to determine the log2foldchange for each overexpression or knockout condition.
Input: Kallisto output directories.
Output: Table with the log2foldchange for each condition.

# Unbiased heatmap
Script name= makeunbiased_heatmap.py

This script takes in the log2foldchange table, merges it with the annotation table, and creates a heatmap of genes that are above or below a log2foldchange cutoff.
Input: log2foldchange table from DESeq2.
Output: PDF with the required heatmap.

# Pathway analysis 
Script name= makeheatmap_NSPpathways.py

This script takes in the log2foldchange table from DESeq2 and a table with a predefined pathway/set of genes in this format: <genename> \t <geneid>.
It creates a heatmap of the desired pathway and outputs a PDF file.
Input: log2foldchange table from DESeq2, predefined set of genes to be plotted.
Output: PDF with the required pathway heatmap.
