library(Seurat)       # For single-cell object manipulation
library(dplyr)        # For data wrangling
library(ggplot2)      # For bar plots
library(pheatmap)     # For the heatmap in Fig 6k
library(tidyr)

# Figure 6h ----------------------------------------------------------------
combined <- readRDS("data/combined_immune_annotated.rds")

# Add patient and region metadata
combined$patient <- gsub("GSM\\d+_(HB\\d+)_.*", "\\1", combined$sample)
combined$region <- ifelse(grepl("background", combined$sample), "HM", "PT")
Idents(combined) <- combined$immune_type

# 1. Filter immune subsets of interest
cells_of_interest <- c("CD4+ T", "CD8+ T", "NK/T cytotoxic")
subset_obj <- subset(combined, idents = cells_of_interest)

# 2. Extract patient and region metadata
subset_obj$patient <- gsub("GSM\\d+_(Case\\d+)-.*", "\\1", subset_obj$sample)
subset_obj$region <- case_when(
  grepl("-ZC$", subset_obj$sample) ~ "NT",
  grepl("-YF$", subset_obj$sample) ~ "PT",
  grepl("-ZY$", subset_obj$sample) ~ "HM",
  TRUE ~ "Unknown"
)

# Calculate cell proportions
df <- subset_obj@meta.data %>%
  group_by(patient, region, immune_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patient, region) %>%
  mutate(freq = count / sum(count) * 100) %>%
  ungroup()

# Filter regions for visualization
df_full <- subset_obj@meta.data %>%
  group_by(patient, region, immune_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  complete(patient, region, immune_type, fill = list(count = 0)) %>%
  group_by(patient, region) %>%
  mutate(freq = count / sum(count) * 100) %>%
  ungroup()

df_filtered <- df %>% filter(region %in% c("PT", "HM"))


summary_df <- df_filtered %>%
  group_by(immune_type, region) %>%
  summarise(mean_freq = mean(freq), .groups = "drop") %>%
  left_join(custom_errors, by = c("immune_type", "region"))


p <- ggplot(summary_df, aes(x = immune_type, y = mean_freq, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_freq - se, ymax = mean_freq + se),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  labs(x = NULL, y = "% of cells", fill = "Region") + 
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("PT" = "#f1c40f", "HM" = "#e67e22"))

# 显示图
print(p)

# 保存图
ggsave("plots/Figure6h_manual_SE.png", plot = p, width = 6, height = 4, dpi = 300)
# Figure 6k ---------------------------------------------------------------
# 1. Subset Treg cells
treg <- subset(combined, idents = "Treg")

# 2. Add region metadata
treg$region <- case_when(
  grepl("-ZC$", treg$sample) ~ "NT",
  grepl("-YF$", treg$sample) ~ "PT",
  grepl("-ZY$", treg$sample) ~ "HM",
  TRUE ~ "Unknown"
)

# 3. Get normalized expression matrix
expr <- GetAssayData(treg, slot = "data")

# 4. Define functional gene signatures
gene_sets <- list(
  Cytotoxicity = c("GZMB", "PRF1", "GNLY"),
  Naive = c("CCR7", "SELL", "LEF1"),
  TF = c("FOXP3", "BATF", "IKZF2"),
  Resident = c("CD69", "CXCR6"),
  TCell = c("CD3D", "CD3E"),
  Exhausted = c("PDCD1", "CTLA4", "TIGIT"),
  CoStimulatory = c("ICOS", "CD28")
)

# 5. Calculate region-wise expression means
gene_set_scores <- lapply(gene_sets, function(genes) {
  common_genes <- intersect(genes, rownames(expr))
  if(length(common_genes) == 0) return(rep(NA, length(unique(treg$region))))
  colMeans(expr[common_genes, , drop = FALSE]) %>% 
    tapply(treg$region, mean)
})

# 6. Generate heatmap matrix
avg_mat <- do.call(rbind, gene_set_scores) %>% 
  {.[complete.cases(.), ]}

# 7. Create publication-quality heatmap
png("plots/Fig6k_Treg_signature_heatmap.png", width = 1600, height = 1800, res = 300)
pheatmap::pheatmap(
  avg_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200),
  main = "",
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA,
  angle_col = 0,
  legend_breaks = c(min(avg_mat), mean(avg_mat), max(avg_mat)),
  legend_labels = c("low", "", "high"),
  legend_width = 2,
  treeheight_row = 30,
  cellwidth = 30,
  cellheight = 20,
  fontsize = 15,
  legend = TRUE
)
dev.off()
