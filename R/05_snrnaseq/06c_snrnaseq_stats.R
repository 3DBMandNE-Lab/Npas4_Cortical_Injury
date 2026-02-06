# 06c_snrnaseq_stats.R
# Calculate snRNA-seq statistics
# Uses the main 81,834 cell snRNA-seq dataset

library(Seurat)
library(dplyr)
source("R/config.R")

# Use the main snRNA-seq dataset
# SNRNASEQ_PATH defined in config.R

cat("Loading snRNA-seq data...\n")
seu <- readRDS(SNRNASEQ_PATH)
cat(sprintf("Loaded: %d cells, %d genes\n", ncol(seu), nrow(seu)))

meta <- seu@meta.data

# Summary by condition
cat("\n=== Cell counts by condition ===\n")
print(table(meta$Condition))

cat("\n=== Cell counts by cell type ===\n")
print(table(meta$celltype_l1))

cat("\n=== Cell type x Condition ===\n")
print(table(meta$celltype_l1, meta$Condition))

# Summary statistics
summary_df <- meta %>%
  group_by(Condition) %>%
  summarize(
    n_cells = n(),
    pct_neurons = 100 * mean(celltype_l1 == "Neuron"),
    pct_microglia = 100 * mean(celltype_l1 == "Microglia"),
    pct_astrocyte = 100 * mean(celltype_l1 == "Astrocyte"),
    pct_oligo = 100 * mean(celltype_l1 == "Oligodendrocyte"),
    .groups = "drop"
  )

cat("\n=== Summary by condition ===\n")
print(summary_df)
write.csv(summary_df, file.path(OUT_TABLES_SNRNASEQ, "snrnaseq_summary.csv"), row.names = FALSE)

# More detailed breakdown
celltype_condition <- meta %>%
  group_by(celltype_l1, Condition) %>%
  summarize(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Condition, values_from = n, values_fill = 0)

cat("\n=== Cells by type and condition ===\n")
print(celltype_condition)
write.csv(celltype_condition, file.path(OUT_TABLES_SNRNASEQ, "snrnaseq_celltype_by_condition.csv"), row.names = FALSE)

cat(sprintf("\nSaved: %s/snrnaseq_*.csv\n", OUT_TABLES_SNRNASEQ))
