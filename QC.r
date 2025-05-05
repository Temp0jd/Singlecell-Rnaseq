library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

# 设置样本数据路径
data_root <- "GSE197177/unzipped_data"
sample_folders <- list.dirs(data_root, full.names = TRUE, recursive = FALSE)

# 初始化列表保存每个 Seurat 对象
seurat_list <- list()

# 逐个样本读取表达矩阵并创建 Seurat 对象
for (folder in sample_folders) {
  sample_name <- basename(folder)
  
  mat <- readMM(file.path(folder, "matrix.mtx"))
  features <- read.delim(file.path(folder, "features.tsv"), header = FALSE)
  barcodes <- read.delim(file.path(folder, "barcodes.tsv"), header = FALSE)
  
  # 设置行名（基因），列名（细胞）
  rownames(mat) <- make.unique(features$V2)  # 用 gene symbol
  colnames(mat) <- barcodes$V1
  
  # 创建 Seurat 对象
  seurat_obj <- CreateSeuratObject(counts = mat, project = sample_name, min.cells = 3, min.features = 200)
  seurat_obj$sample <- sample_name
  seurat_list[[sample_name]] <- seurat_obj
}

# 合并所有样本为一个对象
# 4.1 SCTransform（逐个执行可避免爆内存）
seurat_list <- lapply(seurat_list, function(obj) {
  SCTransform(obj, verbose = FALSE)
})
# 4.2 特征选择 + 整合前准备
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 2000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
# 4.3 寻找锚点
anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT",
                                  anchor.features = features)

# 4.4 执行整合
combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# 4.5 设置 Assay
DefaultAssay(combined) <- "integrated"


-----------------------------------

# 创建 plots 文件夹（如果还不存在）
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Step 2.1: 归一化
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)

# Step 2.2: 识别高变基因（top 2000）
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

# 查看 top10 高变基因（可选）
top10 <- head(VariableFeatures(combined), 10)

# 画出高变基因图
hvf_plot <- VariableFeaturePlot(combined)
hvf_labeled <- LabelPoints(plot = hvf_plot, points = top10, repel = TRUE)

# 保存图像
ggsave("plots/highly_variable_genes.png", plot = hvf_labeled, width = 8, height = 6)

# Step 2.3: 标准化（用于后续 PCA）
combined <- ScaleData(combined, features = VariableFeatures(combined))

# Step 2.4: 主成分分析
combined <- RunPCA(combined, features = VariableFeatures(combined))

# Step 2.5: Elbow Plot 查看主成分拐点
elbow <- ElbowPlot(combined)
ggsave("plots/elbow_plot.png", plot = elbow, width = 5, height = 4)

