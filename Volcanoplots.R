---
  title: "VolcanoPlot"
output: html_notebook
---

# Load necessary library
library(ggrepel)
library(Seurat)
library(tidyverse)
library(ggplot2)

# Assuming `seurat` is your Seurat object
# Here we're grouping by 'orig.ident' to calculate averages for each sample
average_expression <- PseudobulkExpression(seurat, group.by = "orig.ident", assays = "RNA", layer = "data")

# Convert the average expression list to a data frame
average_expression_df <- as.data.frame(average_expression$RNA)
# Add the name "Gene" for the first column
merged_df <- rownames_to_column(average_expression_df , var = "Gene")

# Calculate the average expression for COLD and RT groups
merged_df$avg_COLD <- rowMeans(merged_df[, c("COLD1", "COLD2")])
merged_df$avg_RT <- rowMeans(merged_df[, c("RT1", "RT2")])

# Ensure `orig.ident` has `COLD1`, `COLD2`, `RT1`, and `RT2` as levels
seurat$group <- ifelse(seurat$orig.ident %in% c("COLD1", "COLD2"), "COLD", "RT")

# Set the new grouping as the active identity
Idents(seurat) <- "group"

# Now, run `FindMarkers` comparing the new `COLD` and `RT` groups using DeSeq2
markers <- FindMarkers(seurat, ident.1 = "COLD", ident.2 = "RT", 
                       assay = "RNA", slot = "data", test.use = "DESeq2")

# Add these new columns to our dataframe and add a new column with the -log10 from the pvalue column 

markers <- markers %>%
 mutate(log_pval = -log10(p_val))  
merged_df <- merged_df %>% 
  left_join(markers %>% rownames_to_column(var = "Gene"), by = "Gene")


# Create a column for significance based on thresholds
log2FC_threshold <- 1.5  # threshold for log2 fold change
p_value_threshold <- 0.05  #threshold for p-value

merged_df$Significance <- "Not Significant"
merged_df$Significance[abs(merged_df$avg_log2FC) > log2FC_threshold & merged_df$p_val < p_value_threshold] <- "Significant"

# Specify the gene you want to mark, in case you want to mark one
gene_to_mark <- "gene"

# Find the row index of the gene in merged_df
gene_index <- which(merged_df$Gene == gene_to_mark)

# Create the volcano plot
volcano_plot <- ggplot(merged_df, aes(x = avg_log2FC, y = log_pval)) +
  geom_point(aes(color = Significance), size = 3, alpha = 0.6) +  # Adjust size of circles
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-value") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 12),   # Increase axis text size
    plot.title = element_text(size = 16)   # Increase plot title size
  ) +
  geom_point(data = merged_df[gene_index, ], aes(x = avg_log2FC, y = log_pval), color = "blue", size = 5) + xlim(-4, 4) +  # Set limits for x-axis (adjust as needed)
  ylim(0, 120)  # Set limits for y-axis (adjust as needed)  # Adjust size of specific gene point
  # scale_y_continuous(trans = 'log10') # Set limits for y-axis (adjust as needed)


# Print the plot
print(volcano_plot)

