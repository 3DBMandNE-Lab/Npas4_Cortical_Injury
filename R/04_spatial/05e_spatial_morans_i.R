# 05e_spatial_morans_i.R
# Calculate Moran's I for spatial clustering of signatures
# Input: data/processed/spatial_scored.RDS
# Output: output/tables/spatial/morans_i_statistics.csv

library(dplyr)
source("R/config.R")

# Manual Moran's I calculation (no spdep dependency)
calc_morans_i <- function(x, y, vals, k = 6) {
  n <- length(vals)
  vals <- vals - mean(vals, na.rm = TRUE)

  # Build k-nearest neighbor weights
  coords <- cbind(x, y)
  dists <- as.matrix(dist(coords))
  W <- matrix(0, n, n)

  for (i in 1:n) {
    neighbors <- order(dists[i, ])[2:(k + 1)]  # exclude self
    W[i, neighbors] <- 1 / k
  }

  # Moran's I = (n / S0) * (sum(W * outer(vals, vals)) / sum(vals^2))
  S0 <- sum(W)
  numerator <- sum(W * outer(vals, vals))
  denominator <- sum(vals^2)

  I <- (n / S0) * (numerator / denominator)

  # Expected value under randomness
  E_I <- -1 / (n - 1)

  # Variance (simplified)
  S1 <- 0.5 * sum((W + t(W))^2)
  S2 <- sum((rowSums(W) + colSums(W))^2)
  var_I <- (n * ((n^2 - 3*n + 3) * S1 - n * S2 + 3 * S0^2) -
              sum(vals^4) / (sum(vals^2)^2 / n) * ((n^2 - n) * S1 - 2*n*S2 + 6*S0^2)) /
    ((n - 1) * (n - 2) * (n - 3) * S0^2) - E_I^2

  # Simplify variance calculation
  var_I <- max(var_I, 1e-10)

  z <- (I - E_I) / sqrt(var_I)
  p <- 2 * pnorm(-abs(z))

  list(I = I, E_I = E_I, var = var_I, z = z, p = p)
}

seurat_list <- readRDS(file.path(DATA_PROCESSED, "spatial_scored.RDS"))

results <- list()

for (sid in names(seurat_list)) {
  cat(sprintf("Processing: %s\n", sid))
  meta <- seurat_list[[sid]]@meta.data

  if (!all(c("row", "col") %in% colnames(meta))) {
    cat("  Skipping - no coordinates\n")
    next
  }

  x <- as.numeric(meta$col)
  y <- as.numeric(meta$row)

  # Calculate Moran's I for each signature
  for (score_col in c("Implant_Up", "DAM", "Neuronal", "Complement")) {
    if (!(score_col %in% colnames(meta))) next

    vals <- meta[[score_col]]
    if (all(is.na(vals)) || sd(vals, na.rm = TRUE) == 0) next

    # Replace NA with mean
    vals[is.na(vals)] <- mean(vals, na.rm = TRUE)

    # Moran's I
    mi <- calc_morans_i(x, y, vals, k = 6)

    results[[paste(sid, score_col, sep = "_")]] <- data.frame(
      sample = sid,
      condition = meta$condition[1],
      signature = score_col,
      morans_i = mi$I,
      expected_i = mi$E_I,
      variance = mi$var,
      z_score = mi$z,
      p_value = mi$p
    )

    cat(sprintf("  %s: Moran's I = %.3f (z = %.2f, p = %.2e)\n",
                score_col, mi$I, mi$z, mi$p))
  }
}

results_df <- bind_rows(results)
rownames(results_df) <- NULL
print(results_df)

write.csv(results_df, file.path(OUT_TABLES_SPATIAL, "morans_i_statistics.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %s/morans_i_statistics.csv\n", OUT_TABLES_SPATIAL))
