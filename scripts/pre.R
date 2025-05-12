# Set working directory
setwd("data/GSE197177")  # Replace with your actual path if different

# Create output directory for reformatted data
dir.create("formatted", showWarnings = FALSE)

# List all .gz files in the directory
raw_files <- list.files(pattern = ".gz$", full.names = TRUE)

# Extract sample names by removing suffixes like _barcodes/_features/_matrix
samples <- unique(gsub("_(barcodes|features|matrix).*", "", basename(raw_files)))

# Reorganize files for each sample
for (sample in samples) {
  sample_dir <- file.path("formatted", sample)
  dir.create(sample_dir, showWarnings = FALSE)
  
  sample_files <- raw_files[grepl(sample, raw_files)]
  
  for (f in sample_files) {
    if (grepl("barcodes", f)) {
      file.copy(f, file.path(sample_dir, "barcodes.tsv.gz"))
    }
    if (grepl("features", f)) {
      file.copy(f, file.path(sample_dir, "features.tsv.gz"))
    }
    if (grepl("matrix", f)) {
      file.copy(f, file.path(sample_dir, "matrix.mtx.gz"))
    }
  }
}
