# 05_qc_and_explo_analysis.R
# perform quality control steps and exploratory analysis

# Basic statistics
total_genes <- nrow(counts_raw)
total_samples <- ncol(counts_raw)

# Count statistics per sample
count_stats <- data.frame(
  Sample = colnames(counts_raw),
  Total_Counts = colSums(counts_raw),
  Detected_Genes = colSums(counts_raw > 0),
  stringsAsFactors = FALSE
)
count_stats <- as_tibble(merge(count_stats, metadata, by.x = "Sample", by.y = "row.names"))

#display the count statistics
kable(count_stats, caption = "Count statistics per sample") %>%
  kable_styling(bootstrap_options = "striped")

# display total genes in dataset
cat("Total genes in dataset", total_genes, "\n")

# display total samples in dataset
cat("Total samples:", total_samples, "\n")

# display mean counts per sample
cat("Mean counts per sample:", round(mean(count_stats$Total_Counts), 0), "\n")

# display detected genes per sample
cat("Mean detected genes per sample:", round(mean(count_stats$Detected_Genes), 0), "\n")

# filter out low count genes
keep_genes <- rowSums(counts_raw >= min_count) >= min_samples
counts_filtered <- counts_raw[keep_genes, ]

# export filtered counts matrix as tsv file
write_tsv(counts_filtered, "C:/Users/fraro/Documents/Bioinformatics/chlorella_transcriptomics/data/processed/counts_filtered.tsv")

# display filtering output
cat("Genes before filtering:", nrow(counts_raw), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")
cat("Genes removed:", nrow(counts_raw) - nrow(counts_filtered), "\n")

# create count distribution plots
# first, manipulate data to plot friendly format
count_stats_long <- count_stats %>%
  dplyr::select(Sample, Total_Counts, Detected_Genes) %>%
  pivot_longer(cols = c(Total_Counts, Detected_Genes),
               names_to = "Metric", values_to = "Value")

# create the plot 
p1 <- ggplot(count_stats_long, aes(x = Sample, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Count Statistics per Sample", x = "Sample", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~Metric, scales = "free_y")

# display the plot
print(p1)

# save and export the created plot
ggsave("C:/Users/fraro/Documents/Bioinformatics/chlorella_transcriptomics/output/figures/count_distribution_plots.png", p1)

# Count distribution analysis
# Log2 transformation for visualization
log_counts <- log2(counts_filtered + 1)

# create Density plots, first manipulation in plot friendly form
log_counts_long <- log_counts %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Log2_Count")

# create the plot
p2 <- ggplot(log_counts_long, aes(x = Log2_Count, color = Sample)) +
  geom_density(alpha = 0.7) +
  labs(title = "Distribution of Log2 Counts per Sample",
       x = "Log2(Count + 1)", y = "Density") +
  theme_bw()+
  theme(legend.position = "none")

# display the plot
print(p2)

# save and export the plot
ggsave("C:/Users/fraro/Documents/Bioinformatics/chlorella_transcriptomics/output/figures/distribution_of_log2_counts_per_sample.png", p2)


# create a boxplot of count distributions
# Boxplot of count distributions
p3 <- ggplot(log_counts_long, aes(x = Sample, y = Log2_Count)) +
  geom_boxplot(aes(fill = Sample)) +
  labs(title = "Distribution of Log2 Counts per Sample",
       x = "Sample", y = "Log2(Count + 1)") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# display the plot
print(p3)

# save and export the plot
ggsave("C:/Users/fraro/Documents/Bioinformatics/chlorella_transcriptomics/output/figures/distribution_of_log2_counts_per_sample_boxplot.png", p3)


# create a sample correlation heatmap, correlation after spearman
cor_matrix <- cor(log_counts, method = "spearman")

# Correlation heatmap
heatmap <- pheatmap(cor_matrix,
         main = "Sample Correlation Heatmap (Spearman)",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(0.7, 1, length.out = 101))
print(heatmap)

# save and export the plot
ggsave("C:/Users/fraro/Documents/Bioinformatics/chlorella_transcriptomics/output/figures/sample_corr_heatmap.png", heatmap)


# hierarchical clustering of samples
# Hierarchical clustering of samples
sample_dist <- dist(t(log_counts))
sample_clust <- hclust(sample_dist)

clust_plot <- plot(sample_clust, main = "Sample Clustering",
     xlab = "Samples", ylab = "Distance")

# save and export the plot
ggsave("C:/Users/fraro/Documents/Bioinformatics/chlorella_transcriptomics/output/figures/sample_clustering.png", clust_plot)

# run a principal component analysis
# Remove constant columns (zero variance)
log_counts_novar <- log_counts[apply(log_counts, 1, sd) != 0, ]

# Then run PCA
pca_data <- prcomp(t(log_counts_novar), scale. = TRUE)
pca_summary <- summary(pca_data)

# Calculate variance explained
var_explained <- pca_summary$importance[2, ] * 100

# Create PCA plot data
pca_plot_data <- data.frame(
  Sample = rownames(pca_data$x),
  PC1 = pca_data$x[, 1],
  PC2 = pca_data$x[, 2],
  PC3 = pca_data$x[, 3]
)
pca_plot_data <- merge(pca_plot_data, metadata, by.x = "Sample", by.y = "row.names")

# Color by condition
if ("ammonium" %in% colnames(pca_plot_data)) {
  p4 <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = ammonium)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text_repel(aes(label = Sample)) +
    labs(title = "PCA Plot - Colored by Condition",
         x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
         y = paste0("PC2 (", round(var_explained[2], 1), "%)"))
  
  print(p4)
}

# save and export pca plot
ggsave("C:/Users/fraro/Documents/Bioinformatics/chlorella_transcriptomics/output/figures/pca_condition.png", p4)


# create a scree plot with the variance explained by pca
# Scree plot
scree_data <- data.frame(
  PC = paste0("PC", 1:min(10, ncol(pca_data$x))),
  Variance = var_explained[1:min(10, ncol(pca_data$x))]
)

p5 <- ggplot(scree_data, aes(x = factor(PC, levels = PC), y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "PCA Scree Plot", x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_bw()
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p5)

# save and export the scree plot
ggsave("C:/Users/fraro/Documents/Bioinformatics/chlorella_transcriptomics/output/figures/screeplot_pca_condition.png", p5)
