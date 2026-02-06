# 06b_score_snrnaseq.R
# Add module scores to snRNA-seq data
# Uses the main 81,834 cell snRNA-seq dataset

library(Seurat)
source("R/config.R")

# Use the main snRNA-seq dataset
# SNRNASEQ_PATH defined in config.R

cat("Loading snRNA-seq data...\n")
seu <- readRDS(SNRNASEQ_PATH)
cat(sprintf("Loaded: %d cells, %d genes\n", ncol(seu), nrow(seu)))

# Use RNA assay
DefaultAssay(seu) <- "RNA"

# Define signatures for scoring
signatures <- list(
  Neuronal = c("Rbfox3", "Snap25", "Syt1", "Map2", "Tubb3", "Nefl", "Nefm", "Nefh"),
  Synaptic = c("Gria1", "Gria2", "Grin1", "Shank2", "Homer1", "Dlg4", "Snap25", "Syt1"),
  Activity = c("Npas4", "Arc", "Fos", "Egr1", "Nr4a1", "Jun"),
  Implant_Up = c("C1qa", "C1qb", "C1qc", "Tyrobp", "Trem2", "Gfap", "Vim", "Cd68"),
  DAM = c("Trem2", "Tyrobp", "Apoe", "Lpl", "Cst7", "Cd9", "Ctsd", "Lyz2"),
  Complement = c("C1qa", "C1qb", "C1qc", "C3", "C4b")
)

# Check gene availability and add scores
available <- rownames(seu)
for (sig_name in names(signatures)) {
  genes <- signatures[[sig_name]]
  # Convert to title case for rat genes
  genes_title <- tools::toTitleCase(tolower(genes))
  present <- intersect(genes_title, available)

  cat(sprintf("%s: %d/%d genes available\n", sig_name, length(present), length(genes)))

  if (length(present) >= 3) {
    seu <- AddModuleScore(seu, features = list(present), name = sig_name, seed = 42)
    colnames(seu@meta.data)[ncol(seu@meta.data)] <- sig_name
  }
}

cat("\nCell type distribution:\n")
print(table(seu$celltype_l1))

cat("\nCondition distribution:\n")
print(table(seu$Condition))

cat(sprintf("\n=== snRNA-seq dataset: %d cells, %d genes ===\n", ncol(seu), nrow(seu)))
cat("Module scores added to seu@meta.data\n")
cat("Note: Use readRDS() on the source file directly for analyses\n")
cat(sprintf("Source: %s\n", SNRNASEQ_PATH))
