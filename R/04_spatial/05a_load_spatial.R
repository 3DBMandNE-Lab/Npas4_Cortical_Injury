# 05a_load_spatial.R
# Load Visium data and save as cached RDS
# Output: data/processed/spatial_seurat_list.RDS

library(Seurat)
source("R/config.R")

dir.create("data/processed", showWarnings = FALSE)

samples <- read.csv("data/external/stRNAseq/sample_sheet.csv")
seurat_list <- list()

for (i in 1:nrow(samples)) {
  sid <- samples$sample_id[i]
  cat(sprintf("Loading: %s\n", sid))

  matrix_dir <- file.path("data/external/stRNAseq", sid, "filtered_feature_bc_matrix")
  counts <- Read10X(matrix_dir)
  seu <- CreateSeuratObject(counts, project = sid)
  seu <- NormalizeData(seu, verbose = FALSE)

  # Add spatial coordinates
  pos_file <- file.path("data/external/stRNAseq", sid, "spatial", "tissue_positions.csv")
  if (file.exists(pos_file)) {
    pos <- read.csv(pos_file, header = FALSE)
    colnames(pos) <- c("barcode", "in_tissue", "row", "col", "imagerow", "imagecol")[1:ncol(pos)]
    pos <- pos[pos$barcode %in% colnames(seu), ]
    rownames(pos) <- pos$barcode
    seu$row <- pos[colnames(seu), "row"]
    seu$col <- pos[colnames(seu), "col"]
  }

  # Add metadata
  seu$condition <- samples$condition[i]
  seu$timepoint <- samples$timepoint[i]

  seurat_list[[sid]] <- seu
}

saveRDS(seurat_list, file.path(DATA_PROCESSED, "spatial_seurat_list.RDS"))
cat(sprintf("Saved: %s/spatial_seurat_list.RDS\n", DATA_PROCESSED))
