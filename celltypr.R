library(dplyr)
library(ggplot2)
library(SingleR)
library(celldex)
library(Seurat)

combined <- readRDS("combined_seurat.rds")

# 连接数据层（Seurat v5 新要求）
combined <- JoinLayers(combined)

# 执行差异表达分析
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 保存结果
write.csv(markers, "cluster_markers.csv", row.names = FALSE)

# 提取每个 cluster 的 top5 marker
top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(top5, "cluster_top5_markers.csv", row.names = FALSE)


# 单个基因在UMAP中的表达
FeaturePlot(combined, features = c("CD3D", "MS4A1", "LYZ", "COL1A1", "EPCAM", "PECAM1"), cols = c("lightgrey", "red"))

# DotPlot（推荐）
DotPlot(combined, features = list(
  T_cells = c("CD3D", "IL7R"),
  B_cells = c("MS4A1", "CD79A"),
  Monocytes = c("LYZ", "CD14"),
  NK = c("NKG7", "GNLY"),
  Fibro = c("COL1A1", "PDGFRA"),
  Endo = c("VWF", "PECAM1"),
  Epi = c("EPCAM", "KRT19")
)) + RotatedAxis()


top5 %>%
  group_by(cluster) %>%
  summarise(marker_list = paste(gene, collapse = ", ")) %>%
  print(n = 30)


# 创建注释向量
new.cluster.ids <- c(
  "NK cells", "Enteroendocrine", "Macrophage", "Unknown", "Goblet",
  "Memory T", "CAFs", "Cycling", "Enterocyte", "Enteroendocrine",
  "Epithelial", "Neutrophils", "Macrophage", "Mast cells", "Dying cells",
  "M2 Macrophage", "Endothelial", "Pericyte", "B cells", "Cycling",
  "Cycling", "Plasma cells", "Cytotoxic T", "Fibroblast", "Myofibroblast",
  "Hepatocyte-like"
)

# 应用注释
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
combined$celltype <- Idents(combined)

# 可视化更新

celltype_plot <- DimPlot(combined, group.by = "celltype", label = TRUE, label.size = 4) +
  ggtitle("UMAP with Annotated Cell Types")
ggsave("plots/umap_celltypes_annotated.png", celltype_plot, width = 9, height = 6)



# 切换身份为 cluster
Idents(combined) <- combined$seurat_clusters

# 用 AverageExpression 得到每个 cluster 的平均表达，再做 SingleR 注释
avg_expr <- AverageExpression(combined, return.seurat = FALSE)$RNA
pred_cluster <- SingleR(test = as.matrix(avg_expr), ref = ref, labels = ref$label.fine)

# 应用到 Seurat
cluster_labels <- pred_cluster$labels
names(cluster_labels) <- levels(combined)
combined <- RenameIdents(combined, cluster_labels)
combined$SingleR_cluster <- Idents(combined)

# 绘制 UMAP 图，只显示右侧图例，不在图中加标签
auto_annot_plot <- DimPlot(
  combined,
  group.by = "SingleR_cluster",
  label = FALSE  # 不显示图中文字
) +
  ggtitle("UMAP with SingleR Cluster-Level Annotation")

# 保存图像
ggsave("plots/umap_SingleR_cluster_annotation_legendonly.png", auto_annot_plot, width = 9, height = 6)



