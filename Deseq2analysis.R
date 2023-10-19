'''
Filename: Create_unbiased_heatmap.py
Author: Joel Klein
Date: 2022-07--03
Description: This script creates a heatmap from a deseq2 output file and an annotation file with gene names and
creates a heatmap of the log2foldchange values.
License: CC-BY 4.0 International
Contact: joel.klein@wur.nl
'''

# deseq2 analysis

library("tximport")
library("tximportData")
library("DESeq2")
library("ggplot2")
library("ggrepel")
library("dplyr")
library("svglite")
library("limma")
library("pcaExplorer")


# DEG NSP knockouts
setwd("/path/toworkdir/")

#assigning location of sample information and kallisto output
sample_id <- dir("/path/to/kallistodirectory/kal_out")
sample_id
kall_dirs <- sapply(sample_id, function(id) file.path("/path/to/kallistodirectory/kal_out", id))
kall_dirs

#generate table of sample information
s2c <- read.csv("/path/to/sample_info_2.txt")
s2c <- dplyr::select(s2c, sample, genotype, NP, replication)
s2c <- s2c[order(s2c$sample), ]
s2c <- dplyr::mutate(s2c, path = kall_dirs)
s2c <- s2c[s2c$NP == "NP_Low", ]
s2c

#open and check kallisto files
files <- file.path(s2c$path, "abundance.h5")
files
names(files) <- s2c$sample
files
for (abundance in files) {
  print(paste0(abundance, ifelse(file.exists(abundance), " exists!", " does not exist!")))
}

#Read in the data and create DeSeq2 object
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
head(txi.kallisto$counts)
dds <- DESeqDataSetFromTximport(txi.kallisto,
                                colData = s2c,
                                design= ~ genotype)
dds <- DESeq(dds)

dds$genotype <- relevel(dds$genotype, ref = "C46")

# pcaExplorer(dds)
genotypes <- c('NSP2ox1', 'NSP2ox3', 'NSP2ox6', 'nsp2-9')
res <- list()
for (g in genotypes) {
  res[[g]] <- results(dds, contrast=c('genotype', g, 'C46'), alpha = 0.05)
  print(res[[g]])
}


#generates a table with the log2fold changes NP low
lfc.table <- sapply(res, function(x) x$log2FoldChange)
rownames(lfc.table) <- rownames(res[[1]])
print(lfc.table)
# write the LFC table to file
write.table(cbind(Gene_id = rownames(lfc.table), lfc.table), file = "NSP2_NlowPhigh_LFC.tsv", sep = "\t", row.names = FALSE)



# plot a PCA for our samples
vsd <- vst(dds, blind=FALSE) # perform VST no blind adjustment to size
# rld <- rlog(dds, blind=FALSE) #  rlog transformation
head(assay(vsd), 3)
# plotPCA(vsd, intgroup=c("genotype", "NP"))

# Create PCA plot
pca <- prcomp(t(assay(vsd)))
pca_df <- as.data.frame(pca$x)
pca_df$genotype <- colData(vsd)$genotype
pca_df$NP <- colData(vsd)$NP
pca_df$replicate <- colData(vsd)$replication

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = genotype, fill = genotype)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#F781BF", "#4D4D4D")) +
  scale_shape_manual(values = c(21, 24, 22, 23, 25)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        plot.caption = element_text(size = 12),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  labs(title = "PCA: NP-",
       x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100), "%)"),
       color = "Genotype") +
  geom_text_repel(aes(label = gsub(" ", "", paste(genotype, ".", replicate))), size = 3)

# Remove the legend title
pca_plot <- pca_plot + guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))
# save figure to png
ggsave("pca_plotNlowPhigh.png", pca_plot, width = 8, height = 6, dpi = 300)

