library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

# 1. Load integrated object and set RNA assay
combined <- readRDS("data/combined_log_integrated.rds")
DefaultAssay(combined) <- "RNA"

# 2. Merge all sample data layers (critical step!)
combined[["RNA"]] <- JoinLayers(combined[["RNA"]], layer = "data")

# 3. Standard analysis pipeline
combined <- NormalizeData(combined, assay = "RNA", layer = "data", new.layer = "data")
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, dims = 1:30)

saveRDS(combined, "data/combined_after_clustering.rds")

# 4. Visualize clustering results
DimPlot(combined, label = TRUE, group.by = "seurat_clusters")

# 5. Set cluster IDs as primary identity
Idents(combined) <- "seurat_clusters"

# 6. Identify cluster markers
cluster_markers <- FindAllMarkers(
  combined,
  assay = "RNA",
  only.pos = TRUE,
  logfc.threshold = 0.1,
  min.pct = 0.05
)

# 7. Inspect marker results
dim(cluster_markers)
head(cluster_markers)

# 8. Save marker table (optional)
write.csv(cluster_markers, "data/cluster_markers.csv", row.names = FALSE)

# Extract top 5 markers per cluster
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  select(cluster, gene, avg_log2FC)

# Preview results
print(top_markers, n = 50)

# ---------------------------------------------------------------
# Cell type annotation using SingleR
library(SingleR)
library(celldex)

# Load reference dataset
ref <- HumanPrimaryCellAtlasData()  # Alternative: BlueprintEncodeData()

# Extract expression matrix
expr <- GetAssayData(combined, slot = "data")

# Perform annotation
singler_res <- SingleR(test = expr, ref = ref, labels = ref$label.fine)

# Store results
combined$SingleR <- singler_res$labels

# Map to broader categories
combined$SingleR_group <- case_when(
  grepl("^T\\s?cell", combined$SingleR, ignore.case = TRUE) ~ "T cell",
  grepl("^B\\s?cell", combined$SingleR, ignore.case = TRUE) ~ "B cell",
  grepl("^DC|dendritic", combined$SingleR, ignore.case = TRUE) ~ "Dendritic cell",
  grepl("^Macrophage", combined$SingleR, ignore.case = TRUE) ~ "Macrophage",
  grepl("^Monocyte", combined$SingleR, ignore.case = TRUE) ~ "Monocyte",
  grepl("NK\\s?cell", combined$SingleR, ignore.case = TRUE) ~ "NK cell",
  grepl("^Endothelial", combined$SingleR, ignore.case = TRUE) ~ "Endothelial",
  grepl("^Fibroblast", combined$SingleR, ignore.case = TRUE) ~ "Fibroblast",
  grepl("Stem|iPSC|ESC|embryonic", combined$SingleR, ignore.case = TRUE) ~ "Stem cell",
  grepl("Neutrophil", combined$SingleR, ignore.case = TRUE) ~ "Neutrophil",
  TRUE ~ "Other"
)

# Visualize and save results
p_singler <- DimPlot(combined, 
                     group.by = "SingleR_group", 
                     label = TRUE, 
                     repel = TRUE) +
  ggtitle("UMAP by SingleR Grouped Annotation") +
  theme(legend.text = element_text(size = 8), 
        legend.key.size = unit(0.4, "cm"))

ggsave("plots/umap_by_singler.png", p_singler, width = 12, height = 6, dpi = 300)
saveRDS(singler_res, "data/singler_annotation.rds")

# ---------------------------------------------------------------
# Manual immune cell annotation
# 1. Create cluster-to-celltype mapping
immune_clusters <- c(
  "0" = "CD8+ T", 
  "1" = "CD4+ T", 
  "4" = "DC", 
  "9" = "NK/T cytotoxic",
  "10" = "Treg",
  "11" = "Macrophage",
  "15" = "Neutrophil",
  "20" = "pDC / M2 Macro",
  "21" = "B cell"
)

# 2. Annotate cells
cell_clusters <- as.character(Idents(combined))
immune_type_vec <- immune_clusters[cell_clusters]
immune_type_vec[is.na(immune_type_vec)] <- "Other"
names(immune_type_vec) <- colnames(combined)
combined$immune_type <- immune_type_vec

# Save annotated object
saveRDS(combined, "data/combined_immune_annotated.rds")

# Visualize
p <- DimPlot(combined, 
             group.by = "immune_type", 
             label = TRUE, 
             label.size = 4) + NoLegend()
ggsave("plots/umap_by_immune_type.png", p, width = 8, height = 6, dpi = 300)

# ---------------------------------------------------------------
# Quality control visualizations
# 1. Variable features plot
DefaultAssay(combined) <- "RNA"
VariableFeaturePlot(combined)
ggsave("plots/hvg_variable_feature_plot.png", width = 6, height = 4)

# 2. PCA elbow plot
ElbowPlot(combined, ndims = 50)
ggsave("plots/elbowplot.png", width = 6, height = 4)

# 3. Sample-wise UMAP
DimPlot(combined, group.by = "sample") + ggtitle("By Sample")
ggsave("plots/umap_by_sample.png", width = 6, height = 4)