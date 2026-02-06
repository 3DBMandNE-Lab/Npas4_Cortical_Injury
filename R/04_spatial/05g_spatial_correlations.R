# 05g_spatial_correlations.R
# Pairwise correlations between spatial signatures
# Input: data/processed/spatial_scored.RDS
# Output: output/tables/spatial/pathway_correlations.csv

library(dplyr)
library(tidyr)
source("R/config.R")

seurat_list <- readRDS(file.path(DATA_PROCESSED, "spatial_scored.RDS"))

signatures <- c("Implant_Up", "DAM", "Neuronal", "Complement")

results <- list()

for (sid in names(seurat_list)) {
  meta <- seurat_list[[sid]]@meta.data

  # Get available signatures
  avail <- intersect(signatures, colnames(meta))
  if (length(avail) < 2) next

  # All pairwise correlations
  for (i in 1:(length(avail) - 1)) {
    for (j in (i + 1):length(avail)) {
      sig1 <- avail[i]
      sig2 <- avail[j]

      vals1 <- meta[[sig1]]
      vals2 <- meta[[sig2]]

      valid <- !is.na(vals1) & !is.na(vals2)
      if (sum(valid) < 10) next

      cor_test <- cor.test(vals1[valid], vals2[valid], method = "pearson")

      results[[paste(sid, sig1, sig2, sep = "_")]] <- data.frame(
        sample = sid,
        condition = meta$condition[1],
        signature1 = sig1,
        signature2 = sig2,
        correlation = cor_test$estimate,
        p_value = cor_test$p.value,
        n_spots = sum(valid)
      )
    }
  }
}

results_df <- bind_rows(results)
rownames(results_df) <- NULL

# Pivot to wide format for easier reading
wide_df <- results_df %>%
  mutate(pair = paste(signature1, signature2, sep = "_vs_")) %>%
  select(sample, condition, pair, correlation) %>%
  pivot_wider(names_from = pair, values_from = correlation)

print(results_df)
cat("\n")
print(wide_df)

write.csv(results_df, file.path(OUT_TABLES_SPATIAL, "pathway_correlations.csv"), row.names = FALSE)
write.csv(wide_df, file.path(OUT_TABLES_SPATIAL, "pathway_correlations_wide.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %s/pathway_correlations.csv\n", OUT_TABLES_SPATIAL))
