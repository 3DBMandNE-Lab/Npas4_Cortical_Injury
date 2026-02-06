# 05c_spatial_stats.R
# Calculate spatial statistics (kill zone, correlations)
# Input: data/processed/spatial_scored.RDS
# Output: output/tables/spatial/spatial_stats.csv

library(dplyr)
source("R/config.R")

seurat_list <- readRDS(file.path(DATA_PROCESSED, "spatial_scored.RDS"))

results <- list()

for (sid in names(seurat_list)) {
  meta <- seurat_list[[sid]]@meta.data

  stats <- data.frame(
    sample = sid,
    condition = meta$condition[1],
    n_spots = nrow(meta)
  )

  # Mean scores
  for (col in c("Implant_Up", "DAM", "Neuronal", "Complement")) {
    if (col %in% colnames(meta)) {
      stats[[paste0("mean_", col)]] <- mean(meta[[col]], na.rm = TRUE)
    }
  }

  # Kill zone: high implant + low neuronal
  if (all(c("Implant_Up", "Neuronal") %in% colnames(meta))) {
    thresh_up <- quantile(meta$Implant_Up, 0.75, na.rm = TRUE)
    thresh_neur <- quantile(meta$Neuronal, 0.25, na.rm = TRUE)
    kz <- meta$Implant_Up > thresh_up & meta$Neuronal < thresh_neur
    stats$killzone_pct <- 100 * sum(kz, na.rm = TRUE) / sum(!is.na(kz))
    stats$up_neuronal_cor <- cor(meta$Implant_Up, meta$Neuronal, use = "complete.obs")
  }

  # DAM-Complement correlation
  if (all(c("DAM", "Complement") %in% colnames(meta))) {
    stats$dam_complement_cor <- cor(meta$DAM, meta$Complement, use = "complete.obs")
  }

  results[[sid]] <- stats
}

results_df <- bind_rows(results)
print(results_df)

write.csv(results_df, file.path(OUT_TABLES_SPATIAL, "spatial_stats.csv"), row.names = FALSE)
cat("\nSaved: output/tables/spatial/spatial_stats.csv\n")
