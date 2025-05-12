# 加载所有用到的包
library(CellChat)      # 细胞间通讯分析
library(patchwork)     # 图形排版
library(Seurat)        # 单细胞分析核心工具
library(dplyr)         # 数据操作
library(SingleR)       # 自动细胞注释
library(celldex)       # SingleR 的参考数据集
library(ggplot2)       # 绘图
library(pheatmap)      # 热图绘制
library(Matrix)        # 稀疏矩阵处理
library(stringr)       # 字符串处理
library(mclust)        # 高斯混合模型（用于QC）
library(png)           # PNG图像处理
library(grid)          # 图形系统底层工具

# 保存 sessionInfo() 到文件
sink("sessioninfo.txt")
sessionInfo()
sink()