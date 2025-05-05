# 创建目标目录
dir.create("unzipped_data")

# 获取所有已解压的原始文件
raw_files <- list.files("unzipped_raw", full.names = TRUE)

# 提取样本ID（即 GSM 开头部分）
samples <- unique(gsub("_(matrix|features|barcodes).*", "", basename(raw_files)))

# 对每个样本创建文件夹并整理文件
for (sample in samples) {
  sample_dir <- file.path("unzipped_data", sample)
  dir.create(sample_dir)
  
  # 筛选这个样本的文件
  sample_files <- raw_files[grepl(sample, raw_files)]
  
  # 复制并重命名成 Seurat 所要求的标准名
  for (f in sample_files) {
    if (grepl("matrix", f)) file.copy(f, file.path(sample_dir, "matrix.mtx"))
    if (grepl("features", f)) file.copy(f, file.path(sample_dir, "features.tsv"))
    if (grepl("barcodes", f)) file.copy(f, file.path(sample_dir, "barcodes.tsv"))
  }
}
