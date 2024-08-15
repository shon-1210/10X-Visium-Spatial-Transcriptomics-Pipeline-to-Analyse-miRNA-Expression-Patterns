# Compare multiple clusters together and plot on the same heatmap

# Remember to turn off p-value capping here so we can see how much the p-value has increased.

# Define cluster pairs for FC calculation
# commented out as this is given in the main rmd script to be changed globally.
# cluster_pairs <- list(
#   c("cluster.3", "cluster.6"),
#   c("cluster.4", "cluster.6"),
#   c("cluster.1", "cluster.2")
# )


# Initialize p-value table
pvalTab <- data.frame(row.names = mirs)

# Checking which targets are missing from dataset
sum(!tarTab$`Gene Symbol` %in% rownames(expTab))
numMissingTar <- (length(rownames(expTab)) - (sum(rownames(expTab) %in% tarTab$`Gene Symbol`)))
print(paste("There are", numMissingTar, "missing targets for the chosen miRNAs in the dataset of size", length(rownames(expTab))))

# Keep only those that do have expression
tarTab <- tarTab[tarTab$`Gene Symbol` %in% rownames(expTab), ]

# Keep only those related to mouse (mmu-) miRNAs
tarTab <- tarTab[grepl("^mmu-", tarTab$miRNA), ]

# Order by weighted context++ score, if ties, then context++ score
tarTab <- tarTab[order(tarTab$`weighted context++ score`, tarTab$`context++ score`), ]

# Mean expression values for each cluste already computed earlier
head(meanExpr)

# Loop through each cluster pair
for (specific_clusters in cluster_pairs) {
  
  # Check if specified clusters are in the mean expression data
  if (!all(specific_clusters %in% colnames(meanExpr))) {
    stop("One or more specified clusters not found in the mean expression data.")
  }
  
  # Calculate fold change for the specified clusters
  foldChange <- data.frame(row.names = rownames(meanExpr))
  foldChange[[paste(specific_clusters[1], "vs", specific_clusters[2])]] <- meanExpr[, specific_clusters[1]] - meanExpr[, specific_clusters[2]]
  
  par(mfrow=c(3,2))
  for (mir in mirs) {
    tarGenesAll <- unique(tarTab$`Gene Symbol`[tarTab$miRNA == mir])
    tarGenesTop <- tarGenesAll[seq_len(min(topN, length(tarGenesAll)))]
    
    x <- foldChange[tarGenesTop, 1]
    y <- foldChange[!rownames(foldChange) %in% tarGenesTop, 1]
    pval <- wilcox.test(x, y)$p.value
    pval <- log10(pval)
    pval <- signif(pval, 2) * -sign(median(x) - median(y))
    pvalTab[mir, paste(specific_clusters[1], "vs", specific_clusters[2])] <- signif(pval, 4)
    
    plot(density(y), xlim = c(-2, 2), xlab = "FC density", main = paste(mir, "targets in", paste(specific_clusters[1], "vs", specific_clusters[2])), col = "black")
    lines(density(x), col = 'red')
    legend("topleft", legend = c("non-targets", paste("top", length(x), "targets")), lty = 1, col = c("black", "red"), title = paste("pval =", pval), bty = "n")
  }
}

# Convert the p-value data frame to a matrix
logpvalMatrix <- as.matrix.data.frame(pvalTab)

# Dont Cap values at -10 and 10, to see if it is higher.
# logpvalMatrix[logpvalMatrix > 10] <- 10
# logpvalMatrix[logpvalMatrix < -10] <- -10

# Define the breaks and color palette
breaksList <- seq(-10, 10, by = 0.1)
myColorPalette <- colorRampPalette(c("darkolivegreen3", "white", "red"))(length(breaksList) - 1)

# Generate the heatmap with the transformed values and custom number formatting
pheatmap(logpvalMatrix, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE,
         main = "Heatmap of transformed log10(p-values) for miRNA targets between clusters",
         color = myColorPalette,
         breaks = breaksList,
         fontsize_number = 10,
         angle_col = 45)

