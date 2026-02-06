# 05b_score_spatial.R
# Add module scores to spatial data
# Input: data/processed/spatial_seurat_list.RDS, comparison/implant_specific.csv
# Output: data/processed/spatial_scored.RDS

library(Seurat)
source("R/config.R")

seurat_list <- readRDS(file.path(DATA_PROCESSED, "spatial_seurat_list.RDS"))

# Load signatures
implant <- read.csv(file.path(OUT_TABLES_COMPARISON, "implant_specific.csv"))
genes_up <- implant$gene[implant$category == "Implant_Up"]

SIG_DAM <- c("Spp1", "Tyrobp", "Cd68", "Aif1", "Trem2", "Apoe")
SIG_NEURONAL <- c("Npas4", "Arc", "Fos", "Egr1", "Nr4a1", "Bdnf")
SIG_COMPLEMENT <- c("C1qa", "C1qb", "C1qc", "C3")

for (sid in names(seurat_list)) {
  cat(sprintf("Scoring: %s\n", sid))
  seu <- seurat_list[[sid]]

  for (sig in list(list(genes_up, "Implant_Up"),
                   list(SIG_DAM, "DAM"),
                   list(SIG_NEURONAL, "Neuronal"),
                   list(SIG_COMPLEMENT, "Complement"))) {
    genes <- intersect(sig[[1]], rownames(seu))
    if (length(genes) >= 2) {
      seu <- AddModuleScore(seu, list(genes), name = sig[[2]], nbin = 10)
      colnames(seu@meta.data)[ncol(seu@meta.data)] <- sig[[2]]
    }
  }

  seurat_list[[sid]] <- seu
}

saveRDS(seurat_list, file.path(DATA_PROCESSED, "spatial_scored.RDS"))
cat("Saved: data/processed/spatial_scored.RDS\n")
