# Checking which targets are missing from dataset

# Some gene symbols in miRNA targets might not be in expression
sum(!tarTab$`Gene Symbol` %in% rownames(expTab))

# The number of targets which are missing in the dataset
numMissingTar <- (length(rownames(expTab)) - (sum(rownames(expTab) %in% tarTab$`Gene Symbol`)))
print(paste("There are", numMissingTar, "missing targets for the chosen miRNAs in the dataset of size", length(rownames(expTab))))

dim(tarTab)

# Keep only those that do have expression
tarTab <- tarTab[tarTab$`Gene Symbol` %in% rownames(expTab), ]
dim(tarTab)

# Keep only those related to mouse (mmu-) miRNAs
tarTab <- tarTab[grepl("^mmu-", tarTab$miRNA), ]
dim(tarTab)

# Order by weighted context++ score, if ties, then context++ score
tarTab <- tarTab[order(tarTab$`weighted context++ score`, tarTab$`context++ score`), ]
dim(tarTab)

# Mean expression values for each cluster

# Prepare an empty table to store mean expression values
meanExpr <- data.frame(row.names = rownames(expTab))

for (cluster in sort(unique(cluTab$clusters))) {
  print(paste("Average for cluster:", cluster))
  barcodes <- rownames(cluTab)[cluTab$clusters == cluster]
  # Check if barcodes are valid
  if (length(barcodes) == 0) {
    warning(paste("No barcodes found for cluster:", cluster))
    next
  }
  meanExpr[, cluster] <- rowMeans(expTab[, barcodes, drop = FALSE])
}

# Check if meanExpr is correctly formed and not empty
if (ncol(meanExpr) == 0) {
  stop("meanExpr is empty. Check your clusters and barcodes.")
}

# Ensure that meanExpr contains numeric values
meanExpr <- as.data.frame(lapply(meanExpr, as.numeric), row.names = rownames(meanExpr))

# Check/explore
head(meanExpr)
boxplot(as.matrix(meanExpr), main = "Boxplot of Mean Expressions per Cluster")

# Calculating the FC values for each cluster

# Calculate FC values per cluster
# Create an empty table to store the FC values of each cluster
foldChange <- data.frame(row.names = rownames(meanExpr))

# Looping for each cluster
for (cluster in colnames(meanExpr)) {
  print(cluster)
  foldChange[, cluster] <- meanExpr[, cluster] - rowMeans(meanExpr[, colnames(meanExpr) != cluster, drop = FALSE])
}

# Check/explore
head(foldChange)
boxplot(as.matrix(foldChange), main = "Boxplot of Fold Changes per Cluster")

# Histograms of FC values in each cluster

# Histograms of fold-change per cluster
# change the split according to number of clusters
par(mfrow = mfrow_setting)
for (cluster in colnames(foldChange)) {
  hist(foldChange[, cluster], breaks = 50, main = cluster, xlim = range(-1.5, 1.5))
}

# Plotting density distributions and wilcox test

# Initialize p-value table
pvalTab <- data.frame(row.names = mirs)

# define specific clusters if you want to that should only be present in the final result
#specific_clusters <- c("cluster.3","cluster.4")

# Check miRNA targets for shift in fold-change
# change this according to number of clusters
par(mfrow = mfrow_setting)
for (mir in mirs) {
  print(paste("The miRNA being iterated through now is", mir))
  # We are relying on table already having been sorted above
  tarGenesAll <- unique(tarTab$`Gene Symbol`[tarTab$miRNA == mir])
  # if fewer target genes than topN, choose them all
  tarGenesTop <- tarGenesAll[seq_len(min(topN, length(tarGenesAll)))]
  print(paste("topN value defined earlier was ", topN, " and there are ", length(tarGenesAll), " targets for", mir, "in this dataset."))
  print(paste("The number of targets that will be used for miRNA is", length(tarGenesTop)))
  # Now test these targets for each cluster
  for (cluster in colnames(foldChange)) {       # colnames(foldChange) is for including all clusters # specific_clusters is also an option
    print(cluster)
    x <- foldChange[tarGenesTop, cluster]
    y <- foldChange[!rownames(foldChange) %in% tarGenesTop, cluster] # this needs to be changed for specific cluster comparisons
    targetsMedian <- median(x)
    nonTargetsMedian <- median(y)
    sign_TarvsNonTar <- sign(targetsMedian - nonTargetsMedian)
    print(paste0("The median FC of targets in ", cluster," is ",median(x)))
    print(paste0("The median FC of non-targets in ", cluster," is ",median(y)))
    pval <- wilcox.test(x, y)$p.value
    pval <- log10(pval)
    pval <- signif(pval, 2) * -sign_TarvsNonTar # multiply the sign to know the direction # use - to negate the -ve obtained from log10 transform
    pvalTab[mir, cluster] <- signif(pval, 4) # used to build the heatmap
    plot(density(y), xlim = c(-2, 2), xlab = "FC density", main = paste(mir, "targets in", cluster), col = "black")
    lines(density(x), col = 'red')
    legend("topleft", legend = c("non-targets", paste("top", length(x), "targets")), lty = 1, col = c("black", "red"),
           title = paste("pval =", pval), bty = "n")
  }
}

# Heatmap of p-values and miRNAs in this dataset

# Incorporate the changes Cei said to do
# pvals obtained from wilcox.test need to be log10() transformed before visualising
# change the comments as needed

library(pheatmap)
library(RColorBrewer)


# Convert the p-value data frame to a matrix
logpvalMatrix <- as.matrix.data.frame(pvalTab)

# Cap values at -10 and 10
logpvalMatrix[logpvalMatrix > 10] <- 10
logpvalMatrix[logpvalMatrix < -10] <- -10

# Define the breaks and color palette
breaksList <- seq(-10, 10, by = 0.1)
myColorPalette <- colorRampPalette(c("darkolivegreen3", "white", "red"))(length(breaksList) - 1)

# Here, green p-vals the targets expression is lower than non targets so chances are miRNAs are present in the cluster
# Red p-vals have higher target expression than non-targets so miRNAs are less likely to be present in this cluster

# Generate the heatmap with the transformed values and custom number formatting
pheatmap(logpvalMatrix, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE,
         main = "Heatmap of transformed log10(p-values) for miRNA targets across clusters",
         color = myColorPalette,
         breaks = breaksList,
         fontsize_number = 10,
         angle_col = 45)

