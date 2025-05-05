library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)


# 按照手动注释或 SingleR 注释提取 T cell（你可根据使用的是哪个）
T_cells <- subset(combined, subset = SingleR_cluster %in% c("T_cell:CD4+_effector_memory", "T_cell:gamma-delta", "T_cell:CD8+_effector", "T_cell:Naive"))

#重新标准化和降维（对 T cells 子集）
T_cells <- NormalizeData(T_cells)
T_cells <- FindVariableFeatures(T_cells)
T_cells <- ScaleData(T_cells)
T_cells <- RunPCA(T_cells)

# 可选：看 PCA 拐点
ElbowPlot(T_cells)

#拐点在13左右邻接图 & 聚类（调整 PCA 维度）
T_cells <- FindNeighbors(T_cells, dims = 1:13)
T_cells <- FindClusters(T_cells, resolution = 0.5)  # 可调 0.3~1

#UMAP 可视化（聚类图）
T_cells <- RunUMAP(T_cells, dims = 1:13)

umap_tcell <- DimPlot(T_cells, label = TRUE) +
  ggtitle("T Cell UMAP Clustering (dims=1:13)")

ggsave("plots/umap_Tcell_clustering.png", umap_tcell, width = 8, height = 6)

# 继续差异表达分析 & marker 提取
tcell_markers <- FindAllMarkers(T_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(tcell_markers, "tcell_markers.csv", row.names = FALSE)

tcell_top5 <- tcell_markers %>% group_by(cluster) %>% top_n(5, wt = avg_log2FC)
write.csv(tcell_top5, "tcell_top5_markers.csv", row.names = FALSE)


#T cell marker 可视化
FeaturePlot(T_cells, features = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", "IL7R"))
ggsave("plots/featureplot_Tcells.png", width = 8, height = 6)

DotPlot(T_cells, features = list(
  CD4_T = c("CD4", "IL7R"),
  CD8_T = c("CD8A", "CD8B"),
  Activation = c("GZMB", "PRF1", "IFNG")
)) + RotatedAxis()
ggsave("plots/dotplot_Tcells.png", width = 9, height = 5)



#定义 CD4/CD8/Other 的分组规则
# 添加 CD4/CD8 分群标签
T_cells$T_subtype <- case_when(
  FetchData(T_cells, "CD4")[,1] > 1 & FetchData(T_cells, "CD8A")[,1] < 1 ~ "CD4+ T",
  FetchData(T_cells, "CD8A")[,1] > 1 & FetchData(T_cells, "CD4")[,1] < 1 ~ "CD8+ T",
  FetchData(T_cells, "CD4")[,1] > 1 & FetchData(T_cells, "CD8A")[,1] > 1 ~ "Double Positive",
  TRUE ~ "Other"
)


#统计每个类型的数量和比例
# 统计数量
library(dplyr)
subtype_counts <- T_cells@meta.data %>%
  group_by(T_subtype) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)




#绘制柱状图
barplot <- ggplot(subtype_counts, aes(x = T_subtype, y = percentage, fill = T_subtype)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), vjust = -0.5, size = 4) +
  labs(title = "T Cell Subtype Proportions", x = "Subtype", y = "Percentage") +
  theme_minimal() +
  theme(legend.position = "none")

# 保存图像
ggsave("plots/tcell_subtype_barplot.png", barplot, width = 6, height = 5)
#饼状图


pie_chart <- ggplot(subtype_counts, aes(x = "", y = percentage, fill = T_subtype)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(T_subtype, ": ", percentage, "%")),
            position = position_stack(vjust = 0.5), size = 4) +
  labs(title = "T Cell Subtype Composition") +
  theme_void()

# 保存图像
ggsave("plots/tcell_subtype_piechart.png", pie_chart, width = 6, height = 6)


#------------------------------------------------------------------
#T cell cluster 的 GO 富集分析

top5 <- read.csv("tcell_top5_markers.csv")

top5 %>%
  group_by(cluster) %>%
  summarise(marker_list = paste(gene, collapse = ", ")) %>%
  print(n = 30)



#提取 cluster 5 的 marker genes
cluster5_markers <- FindMarkers(T_cells, ident.1 = 5, min.pct = 0.25, logfc.threshold = 0.25)

# 可选：过滤 logFC 更高的基因
genes <- rownames(cluster5_markers[cluster5_markers$avg_log2FC > 0.5, ])

# 转换为 Entrez ID
library(clusterProfiler)
library(org.Hs.eg.db)

gene_entrez <- bitr(genes,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

# 执行 GO 生物过程（BP）富集分析
ego <- enrichGO(gene = gene_entrez$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",  # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable = TRUE)



# 保存富集分析表格
write.csv(as.data.frame(ego), file = "tcell_cluster5_GO_enrichment.csv", row.names = FALSE)

# 画条形图（前10个通路）
barplot_GO <- barplot(ego, showCategory = 10, title = "GO Enrichment: T cell Cluster 5")
ggsave("plots/tcell_cluster5_GO_barplot.png", barplot_GO, width = 8, height = 6)
