library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)

# Load pre-annotated integrated Seurat object
combined <- readRDS("data/combined_immune_annotated.rds")

combined
head(colnames(combined))
head(combined@meta.data)
table(combined$immune_type)
colnames(combined@meta.data)

# 1. Extract expression matrix (log-normalized)
data.input <- GetAssayData(combined, assay = "RNA", slot = "data")

# 2. Create metadata: use immune_type as grouping variable
meta <- data.frame(group = combined$immune_type, row.names = colnames(combined))

# 3. Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

# 4. Assign ligand–receptor database for human
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# 5. Preprocessing
cellchat <- subsetData(cellchat)  # Filter ligand–receptor pairs
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 6. Compute communication probabilities
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Save CellChat object
saveRDS(cellchat, file = "data/cellchat_recomputed.rds")

# 1. Circle plot of global communication network
png("plots/cellchat_circle.png",
    width  = 800,   # in pixels
    height = 800,
    res    = 150)   # resolution
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale  = TRUE,
                 label.edge    = FALSE)
dev.off()

# 2. Heatmap of signaling network strength
png("plots/cellchat_heatmap.png",
    width  = 1000,
    height = 800,
    res    = 150)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

# 3. Bubble plot of signaling interactions (returns a ggplot object)
bubble_plot <- netVisual_bubble(cellchat,
                                sources.use    = "CD4+ T",
                                targets.use    = NULL,
                                remove.isolate = TRUE,
                                thresh         = 0.1) +
  theme(axis.text.y = element_text(size = 2))

# Save bubble plot
ggsave("plots/cellchat_bubble.png",
       plot   = bubble_plot,
       width  = 6,    # in inches
       height = 8,
       dpi    = 300)
