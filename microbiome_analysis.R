suppressPackageStartupMessages({
  library(vegan)
  library(randomForest)
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
})
dir.create("results", showWarnings = FALSE)

otu_table <- read.csv("data/otu_table.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("data/metadata.csv", stringsAsFactors = FALSE)

prevalence <- colSums(otu_table > 10) / nrow(otu_table)
otu_filtered <- otu_table[, prevalence >= 0.1]
otu_rel <- sweep(otu_filtered, 1, rowSums(otu_filtered), "/")
shannon <- diversity(otu_filtered, index = "shannon")
simpson <- diversity(otu_filtered, index = "simpson")
metadata$Shannon <- shannon[match(metadata$SampleID, rownames(otu_filtered))]
metadata$Simpson <- simpson[match(metadata$SampleID, rownames(otu_filtered))]

shannon_test <- wilcox.test(Shannon ~ Condition, data = metadata)
simpson_test <- wilcox.test(Simpson ~ Condition, data = metadata)

diversity_long <- metadata %>%
  select(SampleID, Condition, Shannon, Simpson) %>%
  pivot_longer(cols = c(Shannon, Simpson), names_to = "Metric", values_to = "Value")

p1 <- ggplot(diversity_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  facet_wrap(~Metric, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Alpha Diversity Comparison",
    y = "Diversity Index", x = ""
    ) + scale_fill_brewer(palette = "Set2") + theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
ggsave("results/alpha_diversity.png", p1, width = 10, height = 5, dpi = 300)
bray_dist <- vegdist(otu_rel, method = "bray")
set.seed(123)
permanova_result <- adonis2(bray_dist ~ Condition, data = metadata, permutations = 999)
pca <- prcomp(otu_rel, scale. = TRUE)
pca_data <- as.data.frame(pca$x[, 1:2])
pca_data$Condition <- metadata$Condition[match(rownames(pca_data), metadata$SampleID)]

# Calculate variance explained
var_explained <- round(100 * summary(pca)$importance[2, 1:2], 1)

p2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(level = 0.95, linetype = 2) +
  theme_bw() +
  labs(
    title = "Beta Diversity (PCA)",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)")
    ) + scale_color_brewer(palette = "Set1") + theme(
  scale_color_brewer(palette = "Set1") +
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
    )
ggsave("results/beta_diversity_pca.png", p2, width = 8, height = 6, dpi = 300)
results_list <- list()

for (taxon in colnames(otu_filtered)) {
  healthy <- otu_filtered[metadata$Condition == "Healthy", taxon]
  disease <- otu_filtered[metadata$Condition == "Cancer", taxon]
  
  mean_healthy <- mean(healthy + 1)
  mean_disease <- mean(disease + 1)
  log2fc <- log2(mean_disease / mean_healthy)
  
  valid_counts <- (sum(healthy) > 0) & (sum(disease) > 0) &
    (length(unique(healthy)) > 1 | length(unique(disease)) > 1)
  
  pvalue <- NA
  if (valid_counts) {
    pvalue <- tryCatch(
      wilcox.test(healthy, disease)$p.value,
      error = function(e) NA
    )
  }
  
  results_list[[taxon]] <- data.frame(
    Taxon = taxon,
    Mean_Healthy = mean(healthy),
    Mean_Cancer = mean(disease),
    Log2FoldChange = log2fc,
    Pvalue = pvalue
  )
}

diff_abundance <- do.call(rbind, results_list)
rownames(diff_abundance) <- NULL

# FDR correction
diff_abundance$FDR <- p.adjust(diff_abundance$Pvalue, method = "BH")

# Sort by p-value
diff_abundance <- diff_abundance %>% arrange(Pvalue)

# Significant taxa
sig_taxa <- diff_abundance %>% filter(FDR < 0.05)
cat(paste("   Found", nrow(sig_taxa), "significantly different taxa (FDR < 0.05)\n"))

# Save results
write.csv(diff_abundance, "results/differential_abundance.csv", row.names = FALSE)

# Plot heatmap of top taxa
if (nrow(sig_taxa) > 0) {
  top_taxa <- head(sig_taxa$Taxon, 15)
  heatmap_data <- otu_rel[, top_taxa, drop = FALSE]
  
  # Annotation
  annotation_row <- data.frame(
    Condition = metadata$Condition,
    row.names = metadata$SampleID
  )
  
  annotation_colors <- list(
    Condition = c(Healthy = "#66C2A5", Cancer = "#FC8D62")
  )
  
  # Scale by column (taxa)
  heatmap_data_scaled <- scale(heatmap_data)
  
  png("results/differential_abundance_heatmap.png", width = 10, height = 8, units = "in", res = 300)
  pheatmap(t(
    heatmap_data_scaled),
    annotation_col = annotation_row,
    annotation_colors = annotation_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    main = "Heatmap of Significantly Different Taxa",
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    fontsize = 10)
  dev.off()
}
rf_predictors <- as.data.frame(otu_rel)
original_taxa_names <- colnames(rf_predictors)
sanitized_taxa_names <- make.names(original_taxa_names, unique = TRUE)
colnames(rf_predictors) <- sanitized_taxa_names

rf_response <- factor(metadata$Condition)

# Train Random Forest model
set.seed(123)
rf_model <- randomForest(
  x = rf_predictors, y = rf_response,
  ntree = 500, importance = TRUE
  )
sink("results/rf_performance.txt")
print(rf_model)
sink()
# Feature importance
importance_mat <- as.data.frame(importance(rf_model))
importance_mat$Sanitized <- rownames(importance_mat)
importance_df <- importance_mat %>%
  mutate(Taxon = original_taxa_names[match(Sanitized, sanitized_taxa_names)]) %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  head(15)

write.csv(importance_df, "results/random_forest_importance.csv", row.names = FALSE)

p3 <- ggplot(importance_df, aes(x = reorder(Taxon, MeanDecreaseAccuracy), 
                                 y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  labs(title = "Top 15 Taxa by Random Forest Importance",
       x = "Taxon", y = "Mean Decrease in Accuracy") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("results/random_forest_importance.png", p3, width = 10, height = 7, dpi = 300)
