#figure复现
#1c 和umap_celltypes_annotated.png类似，但是需要在做修改

unique(combined$sample)

# 构建带 names 的向量：names 是每个细胞的名字
group_vector <- setNames(sampleGroups[combined$sample], colnames(combined))

# 添加为 metadata
combined <- AddMetaData(combined, metadata = group_vector, col.name = "group")





sampleGroups <- c(
  "GSM5910784_Case1-YF" = "PT",
  "GSM5910785_Case1-ZY" = "HM",
  "GSM5910786_Case2-ZC" = "Other",  # ZC 是 normal，不用于分析
  "GSM5910787_Case2-YF" = "PT",
  "GSM5910788_Case2-ZY" = "HM",
  "GSM5910789_Case3-YF" = "PT",
  "GSM5910790_Case3-ZY" = "HM",
  "GSM5910791_Case4-ZY" = "Other"   # 没有配对的 HM
)





# 只保留 PT 和 HM 样本的 T 细胞
T_cells <- subset(combined, subset = group %in% c("PT", "HM") &
                    SingleR_cluster %in% c("T_cell:CD4+_effector_memory", "T_cell:gamma-delta"))

# 统计每种 T cell 亚型在每组中的数量
df <- table(T_cells$group, T_cells$SingleR_cluster) %>% as.data.frame()
colnames(df) <- c("Group", "Subtype", "Count")

# 计算每组的总数，用于归一化
total_counts <- aggregate(Count ~ Group, data = df, sum)
df <- merge(df, total_counts, by = "Group", suffixes = c("", "_Total"))
df$Percentage <- df$Count / df$Count_Total * 100

# 绘图：分组条形图
library(ggpubr)
fig6h <- ggbarplot(df, x = "Subtype", y = "Percentage", fill = "Group",
                   add = "mean_se", position = position_dodge(0.8),
                   palette = c("PT" = "#E7B800", "HM" = "#FC4E07")) +
  rotate() +
  ylab("Proportion (%)") +
  theme_minimal()

print(fig6h)
ggsave("plots/figure6h_Tcell_subtypes_barplot.png", fig6h, width = 7, height = 5)



