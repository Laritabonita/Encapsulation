---
  title: "SpearmanCorrelation"
---
  
# Load necessary libraries
library(Seurat)

# Count the total number of cells per group using metadata
cell_counts <- table(seurat@meta.data$ident)  # 'ident' should be the name of the column that groups the cells
n_cells_COLD <- sum(cell_counts[c("COLD1", "COLD2")], na.rm = TRUE)  # Total cells in COLD
n_cells_RT <- sum(cell_counts[c("RT1", "RT2")], na.rm = TRUE)            # Total cells in RT

# Obtain counts from specific columns
counts <- AggregateExpression(seurat, 
                                 group.by = "ident", 
                                 assays = "RNA", 
                                 slot = "counts", 
                                 return.seurat = FALSE)$RNA  # Assuming "RNA" is the name of the slot

#Calculate the mean of COLD1 and COLD2, and RT1 and RT2
mean_COLD <- rowMeans(counts[, c("COLD1", "COLD2")], na.rm = TRUE)
mean_RT <- rowMeans(counts[, c("RT1", "RT2")], na.rm = TRUE)

# Normalize the mean counts by the total number of cells
mean_COLD_normalized <- mean_COLD / n_cells_COLD
mean_RT_normalized <- mean_RT / n_cells_RT

# Create the plot comparing the normalized means of COLD vs RT
plot(mean_COLD_normalized, mean_RT_normalized,
     xlab = "Mean COLD (Normalized)", ylab = "Mean RT (Normalized)",
     main = paste("Spearman Correlation =", 
                  round(cor(mean_COLD_normalized, mean_RT_normalized, method = "spearman"), 2)),
     xlim = c(0, 1), ylim = c(0, 1),  # Adjust limits as necessary
     pch = 19, col = "black")

# Create the plot comparing the means of COLD vs RT
plot(mean_COLD, mean_RT,
     xlab = "Mean COLD", ylab = "Mean RT",
     main = paste("Spearman Correlation =", 
                  round(cor(mean_COLD_normalized, mean_RT_normalized, method = "spearman"), 2)),
     xlim = c(0, 2000), ylim = c(0, 2000),  # Adjust limits as necessary
     pch = 19, col = "black")

# Add the fitted regression line
abline(lm(mean_RT_normalized ~ mean_COLD_normalized), col = "red", lty = 2)  # Regression line

# Optional: Add a grid for better visualization
grid()

