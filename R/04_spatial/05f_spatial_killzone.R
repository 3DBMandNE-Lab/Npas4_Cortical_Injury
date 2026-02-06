# 05f_spatial_killzone.R
# Kill zone analysis: high inflammation + low neuronal activity
# Input: data/processed/spatial_scored.RDS
# Output: output/tables/spatial/killzone_statistics.csv

library(dplyr)
source("R/config.R")

seurat_list <- readRDS(file.path(DATA_PROCESSED, "spatial_scored.RDS"))

results <- list()

for (sid in names(seurat_list)) {
  meta <- seurat_list[[sid]]@meta.data

  if (!all(c("Implant_Up", "Neuronal") %in% colnames(meta))) {
    cat(sprintf("Skipping %s - missing scores\n", sid))
    next
  }

  # Define kill zone: top 25% Implant_Up AND bottom 25% Neuronal
  thresh_up <- quantile(meta$Implant_Up, 0.75, na.rm = TRUE)
  thresh_neur <- quantile(meta$Neuronal, 0.25, na.rm = TRUE)

  meta$killzone <- meta$Implant_Up > thresh_up & meta$Neuronal < thresh_neur

  # Calculate statistics
  n_spots <- nrow(meta)
  n_killzone <- sum(meta$killzone, na.rm = TRUE)
  pct_killzone <- 100 * n_killzone / n_spots

  # Correlation between Implant_Up and Neuronal
  cor_test <- cor.test(meta$Implant_Up, meta$Neuronal, method = "pearson")

  results[[sid]] <- data.frame(
    sample = sid,
    condition = meta$condition[1],
    n_spots = n_spots,
    n_killzone = n_killzone,
    killzone_pct = pct_killzone,
    up_neuronal_cor = cor_test$estimate,
    up_neuronal_pvalue = cor_test$p.value,
    thresh_up_75 = thresh_up,
    thresh_neur_25 = thresh_neur
  )

  cat(sprintf("%s: Kill zone = %.1f%%, Up-Neuronal r = %.3f\n",
              sid, pct_killzone, cor_test$estimate))
}

results_df <- bind_rows(results)
print(results_df)

write.csv(results_df, file.path(OUT_TABLES_SPATIAL, "killzone_statistics.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %s/killzone_statistics.csv\n", OUT_TABLES_SPATIAL))
