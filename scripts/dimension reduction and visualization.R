library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30)
combined <- RunUMAP(combined, dims = 1:30)
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)
combined$group <- gsub(".*-", "", combined$sample)
p1 <- DimPlot(combined, group.by = "group", label = TRUE) + ggtitle("UMAP by Group")
p2 <- DimPlot(combined, group.by = "seurat_clusters", label = TRUE) + ggtitle("UMAP by Cluster")

ggsave("plots/UMAP_after_integration_group.png", plot = p1, width = 7, height = 6)
ggsave("plots/UMAP_after_integration_cluster.png", plot = p2, width = 7, height = 6)

#before integrated
DefaultAssay(combined) <- "RNA"

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:30)

p3 <- DimPlot(combined, group.by = "sample", label = TRUE) + ggtitle("UMAP Before Integration")
ggsave("plots/UMAP_before_integration.png", plot = p3, width = 7, height = 6)
p4 <- DimPlot(combined, group.by = "seurat_clusters", label = TRUE) + ggtitle("UMAP by Cluster")
ggsave("plots/UMAP_before_integration_cluster.png", plot = p4, width = 7, height = 6)