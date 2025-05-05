# Step 3.1: 构建邻接图
combined <- FindNeighbors(combined, dims = 1:15)

# Step 3.2: 聚类（resolution 可调节聚类数量，建议先 0.5）
combined <- FindClusters(combined, resolution = 0.5)

# Step 3.3: 运行 UMAP 并画图
combined <- RunUMAP(combined, dims = 1:15)

# 保存图像
umap_plot <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 5) +
  ggtitle("UMAP Clustering")

ggsave("plots/umap_clusters.png", plot = umap_plot, width = 8, height = 6)
saveRDS(combined, file = "combined_seurat.rds")
