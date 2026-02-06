# 09a_spatial_density_permutation.R
# Density-controlled permutation test for spatial colocalization
# Controls for UMI density gradients in null model for SPP1-Complement correlation
# Input: data/processed/spatial_scored.RDS or spatial_seurat_list.RDS
# Output: output/tables/spatial/density_permutation_results.csv

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/config.R")

cat("=== Density-Controlled Spatial Permutation Test ===\n\n")

dir.create(file.path(OUT_TABLES_SPATIAL, "permutation"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_FIGURES_SPATIAL, "permutation"), recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load spatial data
# ============================================================

cache_file <- file.path(DATA_PROCESSED, "spatial_seurat_list.RDS")
scored_file <- file.path(DATA_PROCESSED, "spatial_scored.RDS")

if (file.exists(scored_file)) {
  cat("Loading scored spatial data...\n")
  spatial_list <- readRDS(scored_file)
} else if (file.exists(cache_file)) {
  cat("Loading spatial data from cache...\n")
  spatial_list <- readRDS(cache_file)
} else {
  stop("No spatial data found. Run 05a_load_spatial.R first.")
}

# ============================================================
# Define signatures (same as 06k)
# ============================================================

spp1_genes <- c("Spp1", "Gpnmb", "Lgals3", "Fabp5", "Igf1", "Cd63", "Lpl")
complement_genes <- c("C1qa", "C1qb", "C1qc", "C3")

# ============================================================
# Density-controlled permutation function
# ============================================================

#' Permutation test that preserves local density structure
#'
#' @param scores1 Numeric vector of signature 1 scores
#' @param scores2 Numeric vector of signature 2 scores
#' @param umi_counts Numeric vector of UMI counts (density proxy)
#' @param n_perm Number of permutations
#' @param n_bins Number of density bins
#' @return List with observed r, null distribution, and p-value
density_controlled_permutation <- function(scores1, scores2, umi_counts,
                                           n_perm = 1000, n_bins = 10) {

  # Observed correlation
  observed_r <- cor(scores1, scores2, method = "pearson", use = "complete.obs")

  # Create density bins based on UMI counts
  bins <- cut(umi_counts, breaks = n_bins, labels = FALSE)

  # Permutation within bins
  null_r <- numeric(n_perm)

  for (i in seq_len(n_perm)) {
    # Shuffle scores1 WITHIN each density bin
    perm_scores1 <- scores1
    for (b in unique(bins)) {
      idx <- which(bins == b)
      perm_scores1[idx] <- sample(scores1[idx])
    }
    null_r[i] <- cor(perm_scores1, scores2, method = "pearson", use = "complete.obs")
  }

  # Two-sided p-value
  p_value <- mean(abs(null_r) >= abs(observed_r))

  list(
    observed_r = observed_r,
    null_mean = mean(null_r),
    null_sd = sd(null_r),
    null_95ci = quantile(null_r, c(0.025, 0.975)),
    p_value = p_value,
    n_perm = n_perm,
    null_distribution = null_r
  )
}

# ============================================================
# Run permutation test for each sample
# ============================================================

N_PERM <- 1000
set.seed(42)

results_list <- list()
null_distributions <- list()

for (sample_name in names(spatial_list)) {
  cat(sprintf("\n=== Processing %s ===\n", sample_name))
  seu <- spatial_list[[sample_name]]

  # Score if needed
  available_genes <- rownames(seu)

  if (!"SPP1_Module" %in% colnames(seu@meta.data)) {
    spp1_present <- intersect(spp1_genes, available_genes)
    if (length(spp1_present) >= 3) {
      seu <- AddModuleScore(seu, features = list(spp1_present), name = "SPP1_Module", seed = 42)
      colnames(seu@meta.data)[ncol(seu@meta.data)] <- "SPP1_Module"
    }
  }

  if (!"Complement" %in% colnames(seu@meta.data)) {
    comp_present <- intersect(complement_genes, available_genes)
    if (length(comp_present) >= 3) {
      seu <- AddModuleScore(seu, features = list(comp_present), name = "Complement", seed = 42)
      colnames(seu@meta.data)[ncol(seu@meta.data)] <- "Complement"
    }
  }

  # Check if we have the scores
  if (!all(c("SPP1_Module", "Complement") %in% colnames(seu@meta.data))) {
    cat("  Skipping - missing scores\n")
    next
  }

  # Get UMI counts for density proxy
  # Try different possible column names
  if ("nCount_Spatial" %in% colnames(seu@meta.data)) {
    umi_counts <- seu@meta.data$nCount_Spatial
  } else if ("nCount_RNA" %in% colnames(seu@meta.data)) {
    umi_counts <- seu@meta.data$nCount_RNA
  } else {
    # Calculate from counts matrix
    umi_counts <- colSums(GetAssayData(seu, layer = "counts"))
  }

  # Run permutation test
  cat(sprintf("  Running %d permutations...\n", N_PERM))

  perm_result <- density_controlled_permutation(
    scores1 = seu$SPP1_Module,
    scores2 = seu$Complement,
    umi_counts = umi_counts,
    n_perm = N_PERM,
    n_bins = 10
  )

  cat(sprintf("  Observed r = %.3f\n", perm_result$observed_r))
  cat(sprintf("  Null mean  = %.3f (SD = %.3f)\n", perm_result$null_mean, perm_result$null_sd))
  cat(sprintf("  Null 95%% CI = [%.3f, %.3f]\n", perm_result$null_95ci[1], perm_result$null_95ci[2]))
  cat(sprintf("  p-value = %.4f\n", perm_result$p_value))

  # Determine condition from sample name
  condition <- case_when(
    grepl("5B", sample_name) ~ "Control",
    grepl("4A|4B", sample_name) ~ "Acute",
    grepl("7C|8A|8C", sample_name) ~ "Chronic",
    TRUE ~ "Unknown"
  )

  # Store results
  results_list[[sample_name]] <- data.frame(
    sample = sample_name,
    condition = condition,
    n_spots = ncol(seu),
    observed_r = perm_result$observed_r,
    null_mean = perm_result$null_mean,
    null_sd = perm_result$null_sd,
    null_ci_low = perm_result$null_95ci[1],
    null_ci_high = perm_result$null_95ci[2],
    p_value = perm_result$p_value,
    exceeds_null = perm_result$observed_r > perm_result$null_95ci[2]
  )

  null_distributions[[sample_name]] <- data.frame(
    sample = sample_name,
    null_r = perm_result$null_distribution
  )
}

# ============================================================
# Combine and save results
# ============================================================

results_df <- bind_rows(results_list)
null_df <- bind_rows(null_distributions)

cat("\n\n=== SUMMARY: Density-Controlled Permutation Results ===\n")
print(results_df)

write.csv(results_df,
          file.path(OUT_TABLES_SPATIAL, "permutation", "density_permutation_results.csv"),
          row.names = FALSE)

write.csv(null_df,
          file.path(OUT_TABLES_SPATIAL, "permutation", "null_distributions.csv"),
          row.names = FALSE)

# ============================================================
# Visualization
# ============================================================

# Null distribution with observed value for each sample
p_null <- ggplot(null_df, aes(x = null_r)) +
  geom_histogram(bins = 50, fill = "grey70", color = "black", linewidth = 0.2) +
  geom_vline(data = results_df, aes(xintercept = observed_r),
             color = COL_IMPLANT, linewidth = 1, linetype = "dashed") +
  geom_vline(data = results_df, aes(xintercept = null_ci_low),
             color = "grey40", linewidth = 0.5, linetype = "dotted") +
  geom_vline(data = results_df, aes(xintercept = null_ci_high),
             color = "grey40", linewidth = 0.5, linetype = "dotted") +
  facet_wrap(~sample, scales = "free_y", ncol = 2) +
  labs(
    title = "Density-Controlled Permutation Test",
    subtitle = "Red dashed = observed; Grey dotted = 95% null CI",
    x = "SPP1-Complement Correlation (r)",
    y = "Count"
  ) +
  theme_publication()

save_figure(
  file.path(OUT_FIGURES_SPATIAL, "permutation", "null_distributions.png"),
  p_null, width = 10, height = 8
)

# Summary bar plot showing observed vs null range
p_summary <- ggplot(results_df, aes(x = sample)) +
  geom_errorbar(aes(ymin = null_ci_low, ymax = null_ci_high),
                width = 0.3, color = "grey50", linewidth = 0.8) +
  geom_point(aes(y = null_mean), shape = 21, fill = "grey70", size = 3) +
  geom_point(aes(y = observed_r, color = exceeds_null), size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = COL_REF) +
  scale_color_manual(
    values = c("TRUE" = COL_IMPLANT, "FALSE" = COL_NS),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  labs(
    title = "SPP1-Complement Colocalization: Density-Controlled Test",
    subtitle = "Grey = null 95% CI; Colored point = observed",
    x = NULL, y = "Pearson r",
    color = "Exceeds null"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_figure(
  file.path(OUT_FIGURES_SPATIAL, "permutation", "observed_vs_null_summary.png"),
  p_summary, width = 8, height = 5
)

# ============================================================
# Interpretation
# ============================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("DENSITY-CONTROLLED PERMUTATION TEST INTERPRETATION\n")
cat(strrep("=", 70), "\n")

n_sig <- sum(results_df$exceeds_null, na.rm = TRUE)
n_total <- nrow(results_df)

cat(sprintf("\nSamples with r exceeding null 95%% CI: %d / %d\n", n_sig, n_total))

if (n_sig == n_total) {
  cat("\n[RESULT]: CLAIM UPGRADED - Colocalization survives density control\n")
  cat("All samples show SPP1-Complement correlation EXCEEDS what would be\n")
  cat("expected from density alone. The spatial colocalization is NOT\n")
  cat("an artifact of cell density variation.\n")
} else if (n_sig >= n_total / 2) {
  cat("\n[RESULT]: PARTIAL SUPPORT - Majority survive density control\n")
  cat(sprintf("%d/%d samples show real colocalization beyond density.\n", n_sig, n_total))
  cat("Some samples may have density-confounded correlations.\n")
} else {
  cat("\n[RESULT]: CLAIM NOT SUPPORTED - Correlation may be density artifact\n")
  cat("Most correlations do not exceed null expectation after density control.\n")
  cat("The spatial colocalization may reflect density variation, not biology.\n")
}

cat("\n[NEXT STEPS]:\n")
cat("- If significant: Report density-controlled p-values in manuscript\n")
cat("- If not: Revise interpretation of spatial colocalization claim\n")

cat("\n\nOutputs saved to:\n")
cat(sprintf("  %s/permutation/\n", OUT_TABLES_SPATIAL))
cat(sprintf("  %s/permutation/\n", OUT_FIGURES_SPATIAL))
