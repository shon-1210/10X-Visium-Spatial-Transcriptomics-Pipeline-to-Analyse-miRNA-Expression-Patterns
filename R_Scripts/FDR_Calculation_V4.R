# pvalTab is obtained when creating each heatmap
pvalTab[pvalTab > 0] <- -pvalTab[pvalTab > 0]
print(pvalTab)

# Convert back to original pvals
original_pvals <- 10^(pvalTab)
print(original_pvals)

# Convert to vector
numeric_pvaltab <- as.vector(as.matrix(original_pvals))

# Perform the correction, specify n according to the number of p-vals for all tests.
fdrs_with_n <- p.adjust(numeric_pvaltab, method = "BH", n = 672)

print(fdrs_with_n)

# Generate pvalTab back with adjusted pvals
num_cols <- ncol(pvalTab)  # Automatically detect the number of columns
fdrs_with_n_matrix <- matrix(fdrs_with_n, nrow = num_cols, byrow = TRUE)
print(fdrs_with_n_matrix)

# Add row and column names to the matrix
rownames(fdrs_with_n_matrix) <- colnames(original_pvals)
colnames(fdrs_with_n_matrix) <- rownames(original_pvals)
print(fdrs_with_n_matrix)

# Transpose the matrix to match the original_pvals df
fdrs_with_n_matrix_transposed <- t(fdrs_with_n_matrix)

# Print the transposed matrix
print(fdrs_with_n_matrix_transposed)

# Print the original pvals
print(original_pvals)

# Custom function to format numbers for heatmap display
format_numbers <- function(x) {
  ifelse(x < 0.05, signif(x, 4), round(x, 2))
}

# Plot a heatmap of the adjusted pvals
# Define the breaks and color palette
# Ensure unique breaks by using a small increment
breaksList <- c(seq(0.000001, 0.05, by = 0.001))
myColorPalette <- colorRampPalette(c("red", "lightblue", "white"))(length(breaksList) - 1)

# Generate the heatmap with the transformed values and custom number formatting
pheatmap(fdrs_with_n_matrix_transposed, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = matrix(format_numbers(fdrs_with_n_matrix_transposed), 
                                  nrow = nrow(fdrs_with_n_matrix_transposed), 
                                  ncol = ncol(fdrs_with_n_matrix_transposed)),
         main = "Heatmap of adjusted p values for miRNA targets between clusters",
         color = myColorPalette,
         breaks = breaksList,
         fontsize_number = 10,
         angle_col = 45)
