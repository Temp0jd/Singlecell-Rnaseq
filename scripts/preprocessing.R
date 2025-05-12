library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

# Set sample data paths
data_root <- "data/GSE197177/formatted"
sample_folders <- list.dirs(data_root, full.names = TRUE, recursive = FALSE)

# Initialize list to store Seurat objects
seurat_list <- list()

# Read expression matrices and create Seurat objects for each sample
for (folder in sample_folders) {
  sample_name <- basename(folder)
  
  # Load matrix files
  mat <- readMM(file.path(folder, "matrix.mtx.gz"))
  features <- read.delim(file.path(folder, "features.tsv.gz"), header = FALSE)
  barcodes <- read.delim(file.path(folder, "barcodes.tsv.gz"), header = FALSE)
  
  # Set row names (genes) and column names (cells)
  rownames(mat) <- make.unique(features$V2)  # Use gene symbols
  colnames(mat) <- barcodes$V1
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = mat,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  seurat_obj$sample <- sample_name
  seurat_list[[sample_name]] <- seurat_obj
}

# Preprocessing pipeline (Normalization, Variable Features, Scaling, PCA)
seurat_list <- lapply(seurat_list, function(obj) {
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 1500)
  obj <- ScaleData(obj, features = VariableFeatures(obj))
  obj <- RunPCA(obj, features = VariableFeatures(obj))
  return(obj)
})

saveRDS(seurat_list, file = "data/seurat_list_normalized.rds")

# Integration workflow
# 1. Select integration features
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 1500)
saveRDS(features, file = "data/integration_features.rds")

# 2. Find integration anchors
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  anchor.features = features,
  reduction = "rpca"
)
saveRDS(anchors, file = "data/integration_anchors.rds")

# 3. Perform integration
combined <- IntegrateData(anchorset = anchors)
DefaultAssay(combined) <- "integrated"

# 4. Save final integrated object
saveRDS(combined, file = "data/combined_log_integrated.rds")

# ---------------------------
# Quality Control Section
library(stringr)
library(mclust)
dir.create("plots/qc", showWarnings = FALSE)

# Initialize QC summary table
qc_summary <- data.frame(
  Sample = character(),
  Cells_Before = numeric(),
  Cells_After = numeric(),
  Genes_Before = numeric(),
  Genes_After = numeric(),
  nFeature_Gaussians = numeric(),
  stringsAsFactors = FALSE
)

for (sample in names(seurat_list)) {
  obj <- seurat_list[[sample]]
  
  # Create short sample labels
  short_sample <- str_replace(sample, "GSM\\d+_", "")
  obj$short_label <- short_sample
  Idents(obj) <- "short_label"
  
  # Calculate mitochondrial percentage
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  seurat_list[[sample]] <- obj
  
  # Generate QC violin plots
  p1 <- VlnPlot(
    obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0.01
  ) +
    plot_annotation(title = paste("QC Metrics -", short_sample)) +
    theme(
      axis.text.x = element_text(size = 8),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 12)
    )
  
  ggsave(
    filename = paste0("plots/qc/QC_", short_sample, ".png"),
    plot = p1,
    width = 15,
    height = 4,
    dpi = 300
  )
  
  # Collect QC metrics
  before_cells <- ncol(obj)
  before_genes <- nrow(obj)
  
  # Apply basic filters
  filtered <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  after_cells <- ncol(filtered)
  after_genes <- nrow(filtered)
  
  # Bimodal detection: Fit GMM to nFeature distribution
  x <- obj$nFeature_RNA
  x <- x[x > 0]
  gmm <- Mclust(x, G = 1:3)
  best_G <- gmm$G
  
  # Save GMM visualization
  png(
    paste0("plots/qc/GMM_nFeature_", short_sample, ".png"),
    width = 800,
    height = 600
  )
  hist(x, breaks = 100, freq = FALSE, main = paste("GMM fit -", short_sample), xlab = "nFeature_RNA")
  lines(density(x), col = "black")
  for (k in 1:best_G) {
    lines(density(gmm$data[gmm$classification == k]), col = k + 1, lwd = 2)
  }
  dev.off()
  
  # Update summary table
  qc_summary <- rbind(qc_summary, data.frame(
    Sample = short_sample,
    Cells_Before = before_cells,
    Cells_After = after_cells,
    Genes_Before = before_genes,
    Genes_After = after_genes,
    nFeature_Gaussians = best_G
  ))
}

# Generate combined QC plots
library(png)
library(grid)
library(patchwork)

# Combine GMM plots
gmm_files <- list.files("plots/qc", pattern = "^GMM_nFeature_.*\\.png$", full.names = TRUE)

plot_list <- lapply(gmm_files, function(file) {
  img <- readPNG(file)
  grob <- rasterGrob(img, interpolate = TRUE)
  ggplot() + annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    theme_void() + ggtitle(basename(file))
})

combined_gmm_plot <- wrap_plots(plotlist = plot_list, ncol = 3)
ggsave(
  "plots/qc/ALL_GMM_combined.png",
  combined_gmm_plot,
  width = 30,
  height = ceiling(length(plot_list) / 3) * 5,
  dpi = 300
)

# Generate combined violin plots
qc_plots <- lapply(seurat_list, function(obj) {
  short_sample <- unique(obj$short_label)
  VlnPlot(
    obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0.02
  ) +
    plot_annotation(title = short_sample) +
    theme(
      axis.text.x = element_text(size = 7),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 10)
    )
})

wrap_plot <- wrap_plots(qc_plots, ncol = 2)
ggsave(
  "plots/qc/QC_all_samples_combined.png",
  wrap_plot,
  width = 20,
  height = ceiling(length(qc_plots) / 2) * 4,
  dpi = 600
)

# Save QC summary
write.csv(qc_summary, "data/qc_summary.csv", row.names = FALSE)