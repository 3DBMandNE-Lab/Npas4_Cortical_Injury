# 06a_load_snrnaseq.R
# Load full snRNA-seq dataset
# Input: Seurat object (path set by SNRNASEQ_PATH in config.R)
# Output: Seurat object in memory, metadata summary

library(Seurat)
source("R/config.R")

SNRNASEQ_SOURCE <- SNRNASEQ_PATH  # defined in config.R

cat("Loading full snRNA-seq dataset...\n")
seu <- readRDS(SNRNASEQ_SOURCE)

cat(sprintf("Loaded: %d cells, %d genes\n", ncol(seu), nrow(seu)))

cat("\nMetadata columns:\n")
print(colnames(seu@meta.data))

cat("\nConditions:\n")
if ("Condition" %in% colnames(seu@meta.data)) {
  print(table(seu$Condition))
}

cat("\nCell types (L1):\n")
if ("celltype_l1" %in% colnames(seu@meta.data)) {
  print(table(seu$celltype_l1))
}

cat("\nCell types (L2):\n")
if ("celltype_l2" %in% colnames(seu@meta.data)) {
  print(table(seu$celltype_l2))
}

# Save metadata for inspection
dir.create(DATA_PROCESSED, showWarnings = FALSE, recursive = TRUE)
write.csv(seu@meta.data, file.path(DATA_PROCESSED, "snrnaseq_meta.csv"))
cat(sprintf("\nSaved metadata: %s/snrnaseq_meta.csv\n", DATA_PROCESSED))
