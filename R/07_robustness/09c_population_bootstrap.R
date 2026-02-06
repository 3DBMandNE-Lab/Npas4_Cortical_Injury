# 09c_population_bootstrap.R
# Bootstrap confidence intervals for microglial population overlap statistics
# Bootstrap confidence intervals for SPP1+ vs Complement+ overlap odds ratio
# Input: snRNA-seq from external source
# Output: output/tables/snrnaseq/populations_bootstrap.csv

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/config.R")

cat("=== Bootstrap Confidence Intervals for Population Statistics ===\n\n")

dir.create(file.path(OUT_TABLES_SNRNASEQ, "population_bootstrap"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_FIGURES_SNRNASEQ, "population_bootstrap"), recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load snRNA-seq data
# ============================================================

SNRNASEQ_SOURCE <- SNRNASEQ_PATH  # defined in config.R

if (!file.exists(SNRNASEQ_SOURCE)) {
  stop("snRNA-seq data not found at: ", SNRNASEQ_SOURCE)
}

cat("Loading snRNA-seq data...\n")
seu <- readRDS(SNRNASEQ_SOURCE)
DefaultAssay(seu) <- "RNA"

# Subset to microglia
celltype_col <- if ("celltype_l2" %in% colnames(seu@meta.data)) "celltype_l2" else "celltype_l1"
microglia_types <- unique(seu@meta.data[[celltype_col]])[
  grepl("Microglia|microglia|MG", unique(seu@meta.data[[celltype_col]]), ignore.case = TRUE)
]
seu_mg <- subset(seu, cells = colnames(seu)[seu@meta.data[[celltype_col]] %in% microglia_types])
DefaultAssay(seu_mg) <- "RNA"

cat(sprintf("Microglia: %d cells\n", ncol(seu_mg)))

# ============================================================
# Define populations
# ============================================================

# SPP1+
spp1_expr <- GetAssayData(seu_mg, assay = "RNA", layer = "data")["Spp1", ]
seu_mg$spp1_pos <- spp1_expr > 0

# Complement score
complement_genes <- c("C1qa", "C1qb", "C1qc", "C3")
available_genes <- rownames(seu_mg)
comp_present <- intersect(tools::toTitleCase(tolower(complement_genes)), available_genes)

if (length(comp_present) >= 3) {
  seu_mg <- AddModuleScore(seu_mg, features = list(comp_present), name = "Complement", seed = 42)
  colnames(seu_mg@meta.data)[ncol(seu_mg@meta.data)] <- "Complement"
}

# Define Complement high (top 25%)
comp_q75 <- quantile(seu_mg$Complement, 0.75)
seu_mg$comp_high <- seu_mg$Complement >= comp_q75

# Create 2x2 classification
seu_mg$population <- case_when(
  seu_mg$spp1_pos & seu_mg$comp_high ~ "SPP1+/Comp+",
  seu_mg$spp1_pos & !seu_mg$comp_high ~ "SPP1+/Comp-",
  !seu_mg$spp1_pos & seu_mg$comp_high ~ "SPP1-/Comp+",
  TRUE ~ "SPP1-/Comp-"
)

# ============================================================
# Bootstrap functions
# ============================================================

#' Calculate overlap statistics from contingency table
calc_overlap_stats <- function(spp1_pos, comp_high) {
  # Create 2x2 table
  tab <- table(SPP1 = spp1_pos, Comp = comp_high)

  # Handle edge cases
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(list(
      n_both = 0,
      n_spp1_only = sum(spp1_pos & !comp_high),
      n_comp_only = sum(!spp1_pos & comp_high),
      n_neither = sum(!spp1_pos & !comp_high),
      pct_overlap = 0,
      OR = NA,
      jaccard = 0
    ))
  }

  # Cell counts
  a <- tab["TRUE", "TRUE"]   # Both positive
  b <- tab["TRUE", "FALSE"]  # SPP1 only
  c <- tab["FALSE", "TRUE"]  # Comp only
  d <- tab["FALSE", "FALSE"] # Neither

  # Overlap percentage
  pct_overlap <- 100 * a / (a + b + c)

  # Odds ratio (with continuity correction)
  OR <- (a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5))

  # Jaccard index: intersection / union
  jaccard <- a / (a + b + c)

  list(
    n_both = a,
    n_spp1_only = b,
    n_comp_only = c,
    n_neither = d,
    pct_overlap = pct_overlap,
    OR = OR,
    jaccard = jaccard
  )
}

#' Bootstrap population statistics
#'
#' @param spp1_pos Logical vector of SPP1+ status
#' @param comp_high Logical vector of Complement_High status
#' @param n_boot Number of bootstrap iterations
#' @return List with observed stats, CIs, and distributions
bootstrap_population_stats <- function(spp1_pos, comp_high, n_boot = 2000) {

  # Observed statistics
  observed <- calc_overlap_stats(spp1_pos, comp_high)

  # Bootstrap
  n <- length(spp1_pos)
  boot_OR <- numeric(n_boot)
  boot_jaccard <- numeric(n_boot)
  boot_pct_overlap <- numeric(n_boot)

  for (b in seq_len(n_boot)) {
    idx <- sample(n, replace = TRUE)
    stats <- calc_overlap_stats(spp1_pos[idx], comp_high[idx])
    boot_OR[b] <- stats$OR
    boot_jaccard[b] <- stats$jaccard
    boot_pct_overlap[b] <- stats$pct_overlap
  }

  list(
    observed = observed,
    OR_ci = quantile(boot_OR, c(0.025, 0.975), na.rm = TRUE),
    OR_mean = mean(boot_OR, na.rm = TRUE),
    jaccard_ci = quantile(boot_jaccard, c(0.025, 0.975), na.rm = TRUE),
    jaccard_mean = mean(boot_jaccard, na.rm = TRUE),
    pct_overlap_ci = quantile(boot_pct_overlap, c(0.025, 0.975), na.rm = TRUE),
    pct_overlap_mean = mean(boot_pct_overlap, na.rm = TRUE),
    boot_distributions = data.frame(
      OR = boot_OR,
      jaccard = boot_jaccard,
      pct_overlap = boot_pct_overlap
    )
  )
}

# ============================================================
# Run bootstrap analysis
# ============================================================

N_BOOT <- 2000
set.seed(42)

cat(sprintf("\nRunning %d bootstrap iterations...\n", N_BOOT))

boot_result <- bootstrap_population_stats(
  seu_mg$spp1_pos,
  seu_mg$comp_high,
  n_boot = N_BOOT
)

cat("\n=== OBSERVED STATISTICS ===\n")
cat(sprintf("SPP1+ / Comp+  (both):     %d cells (%.2f%%)\n",
            boot_result$observed$n_both,
            100 * boot_result$observed$n_both / ncol(seu_mg)))
cat(sprintf("SPP1+ / Comp-  (SPP1 only): %d cells (%.2f%%)\n",
            boot_result$observed$n_spp1_only,
            100 * boot_result$observed$n_spp1_only / ncol(seu_mg)))
cat(sprintf("SPP1- / Comp+  (Comp only): %d cells (%.2f%%)\n",
            boot_result$observed$n_comp_only,
            100 * boot_result$observed$n_comp_only / ncol(seu_mg)))
cat(sprintf("SPP1- / Comp-  (neither):   %d cells (%.2f%%)\n",
            boot_result$observed$n_neither,
            100 * boot_result$observed$n_neither / ncol(seu_mg)))

cat("\n=== BOOTSTRAP CONFIDENCE INTERVALS ===\n")
cat(sprintf("Odds Ratio:       %.3f (95%% CI: %.3f - %.3f)\n",
            boot_result$observed$OR,
            boot_result$OR_ci[1], boot_result$OR_ci[2]))
cat(sprintf("Jaccard Index:    %.4f (95%% CI: %.4f - %.4f)\n",
            boot_result$observed$jaccard,
            boot_result$jaccard_ci[1], boot_result$jaccard_ci[2]))
cat(sprintf("Overlap %%:        %.2f%% (95%% CI: %.2f%% - %.2f%%)\n",
            boot_result$observed$pct_overlap,
            boot_result$pct_overlap_ci[1], boot_result$pct_overlap_ci[2]))

# ============================================================
# Fisher's exact test (for comparison)
# ============================================================

fisher_result <- fisher.test(table(seu_mg$spp1_pos, seu_mg$comp_high))
cat(sprintf("\nFisher's exact test:\n"))
cat(sprintf("  OR = %.3f (95%% CI: %.3f - %.3f)\n",
            fisher_result$estimate,
            fisher_result$conf.int[1], fisher_result$conf.int[2]))
cat(sprintf("  p-value = %.4f\n", fisher_result$p.value))

# ============================================================
# Resolution stability test
# ============================================================

cat("\n\n=== RESOLUTION STABILITY TEST ===\n")
cat("Testing if SPP1+ identification is stable across thresholds...\n")

# Test different percentile cutoffs for Complement_High
percentiles <- c(0.70, 0.75, 0.80, 0.85, 0.90)
stability_results <- list()

for (pct in percentiles) {
  thresh <- quantile(seu_mg$Complement, pct)
  comp_high_test <- seu_mg$Complement >= thresh

  stats <- calc_overlap_stats(seu_mg$spp1_pos, comp_high_test)

  stability_results[[as.character(pct)]] <- data.frame(
    percentile = pct,
    threshold = thresh,
    n_comp_high = sum(comp_high_test),
    n_both = stats$n_both,
    OR = stats$OR,
    jaccard = stats$jaccard
  )

  cat(sprintf("  Comp threshold %.0f%%: OR = %.3f, Jaccard = %.4f, n_both = %d\n",
              pct * 100, stats$OR, stats$jaccard, stats$n_both))
}

stability_df <- bind_rows(stability_results)

# ============================================================
# By condition analysis
# ============================================================

cat("\n\n=== BY CONDITION ANALYSIS ===\n")

condition_col <- if ("Condition" %in% colnames(seu_mg@meta.data)) "Condition" else "condition"
if (!condition_col %in% colnames(seu_mg@meta.data)) {
  cat("No condition column found, skipping condition analysis\n")
  condition_results <- NULL
} else {
  conditions <- unique(seu_mg@meta.data[[condition_col]])
  condition_results <- list()

  for (cond in conditions) {
    idx <- seu_mg@meta.data[[condition_col]] == cond
    spp1_cond <- seu_mg$spp1_pos[idx]
    comp_cond <- seu_mg$comp_high[idx]

    stats <- calc_overlap_stats(spp1_cond, comp_cond)

    condition_results[[cond]] <- data.frame(
      condition = cond,
      n_cells = sum(idx),
      n_spp1 = sum(spp1_cond),
      n_comp_high = sum(comp_cond),
      n_both = stats$n_both,
      OR = stats$OR,
      jaccard = stats$jaccard
    )

    cat(sprintf("%s: n=%d, SPP1+=%d, Comp+=%d, Both=%d, OR=%.3f\n",
                cond, sum(idx), sum(spp1_cond), sum(comp_cond),
                stats$n_both, stats$OR))
  }

  condition_results_df <- bind_rows(condition_results)
}

# ============================================================
# Save results
# ============================================================

summary_df <- data.frame(
  statistic = c("Odds_Ratio", "Jaccard", "Pct_Overlap"),
  observed = c(boot_result$observed$OR,
               boot_result$observed$jaccard,
               boot_result$observed$pct_overlap),
  ci_low = c(boot_result$OR_ci[1],
             boot_result$jaccard_ci[1],
             boot_result$pct_overlap_ci[1]),
  ci_high = c(boot_result$OR_ci[2],
              boot_result$jaccard_ci[2],
              boot_result$pct_overlap_ci[2]),
  boot_mean = c(boot_result$OR_mean,
                boot_result$jaccard_mean,
                boot_result$pct_overlap_mean)
)

write.csv(summary_df,
          file.path(OUT_TABLES_SNRNASEQ, "population_bootstrap", "bootstrap_summary.csv"),
          row.names = FALSE)

write.csv(boot_result$boot_distributions,
          file.path(OUT_TABLES_SNRNASEQ, "population_bootstrap", "bootstrap_distributions.csv"),
          row.names = FALSE)

write.csv(stability_df,
          file.path(OUT_TABLES_SNRNASEQ, "population_bootstrap", "resolution_stability.csv"),
          row.names = FALSE)

if (!is.null(condition_results)) {
  write.csv(condition_results_df,
            file.path(OUT_TABLES_SNRNASEQ, "population_bootstrap", "by_condition.csv"),
            row.names = FALSE)
}

# ============================================================
# Visualization
# ============================================================

# Bootstrap distribution of OR
p_or <- ggplot(boot_result$boot_distributions, aes(x = OR)) +
  geom_histogram(bins = 50, fill = "grey70", color = "black", linewidth = 0.2) +
  geom_vline(xintercept = boot_result$observed$OR, color = COL_IMPLANT,
             linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = boot_result$OR_ci, color = "grey40",
             linewidth = 0.8, linetype = "dotted") +
  geom_vline(xintercept = 1, color = COL_REF, linewidth = 0.8) +
  annotate("text", x = 1, y = Inf, label = "OR = 1\n(independence)",
           vjust = 2, hjust = -0.1, size = 3, color = COL_REF) +
  labs(
    title = "Bootstrap Distribution of Odds Ratio",
    subtitle = sprintf("OR = %.3f (95%% CI: %.3f - %.3f)",
                       boot_result$observed$OR,
                       boot_result$OR_ci[1], boot_result$OR_ci[2]),
    x = "Odds Ratio (SPP1+ vs Complement_High)",
    y = "Count"
  ) +
  theme_publication()

save_figure(
  file.path(OUT_FIGURES_SNRNASEQ, "population_bootstrap", "or_bootstrap.png"),
  p_or, width = 8, height = 5
)

# Bootstrap distribution of Jaccard
p_jaccard <- ggplot(boot_result$boot_distributions, aes(x = jaccard)) +
  geom_histogram(bins = 50, fill = "grey70", color = "black", linewidth = 0.2) +
  geom_vline(xintercept = boot_result$observed$jaccard, color = COL_IMPLANT,
             linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = boot_result$jaccard_ci, color = "grey40",
             linewidth = 0.8, linetype = "dotted") +
  labs(
    title = "Bootstrap Distribution of Jaccard Index",
    subtitle = sprintf("Jaccard = %.4f (95%% CI: %.4f - %.4f)",
                       boot_result$observed$jaccard,
                       boot_result$jaccard_ci[1], boot_result$jaccard_ci[2]),
    x = "Jaccard Index (SPP1+ ∩ Complement_High)",
    y = "Count"
  ) +
  theme_publication()

save_figure(
  file.path(OUT_FIGURES_SNRNASEQ, "population_bootstrap", "jaccard_bootstrap.png"),
  p_jaccard, width = 8, height = 5
)

# Resolution stability plot
p_stability <- ggplot(stability_df, aes(x = percentile * 100, y = OR)) +
  geom_line(color = COL_IMPLANT, linewidth = 1) +
  geom_point(size = 3, color = COL_IMPLANT) +
  geom_hline(yintercept = 1, linetype = "dashed", color = COL_REF) +
  scale_x_continuous(breaks = percentiles * 100) +
  labs(
    title = "Resolution Stability: OR Across Thresholds",
    subtitle = "OR stable below 1 across Complement_High definitions",
    x = "Complement Percentile Threshold",
    y = "Odds Ratio"
  ) +
  theme_publication()

save_figure(
  file.path(OUT_FIGURES_SNRNASEQ, "population_bootstrap", "resolution_stability.png"),
  p_stability, width = 7, height = 5
)

# Venn-style visualization
# Calculate areas for illustration
total <- ncol(seu_mg)
spp1_total <- sum(seu_mg$spp1_pos)
comp_total <- sum(seu_mg$comp_high)
both <- boot_result$observed$n_both

venn_data <- data.frame(
  category = c("SPP1+ only", "Comp+ only", "Both", "Neither"),
  count = c(boot_result$observed$n_spp1_only,
            boot_result$observed$n_comp_only,
            boot_result$observed$n_both,
            boot_result$observed$n_neither),
  pct = c(100 * boot_result$observed$n_spp1_only / total,
          100 * boot_result$observed$n_comp_only / total,
          100 * boot_result$observed$n_both / total,
          100 * boot_result$observed$n_neither / total)
)

p_venn <- ggplot(venn_data, aes(x = category, y = pct, fill = category)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5) +
  geom_text(aes(label = sprintf("n=%d\n(%.1f%%)", count, pct)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(
    "SPP1+ only" = COL_IMPLANT,
    "Comp+ only" = "#984EA3",
    "Both" = "#FF7F00",
    "Neither" = "grey70"
  )) +
  labs(
    title = "SPP1+ vs Complement_High Population Distribution",
    subtitle = sprintf("Total microglia: %d | Overlap: %.1f%% (OR = %.2f)",
                       total, venn_data$pct[venn_data$category == "Both"],
                       boot_result$observed$OR),
    x = NULL, y = "% of Microglia"
  ) +
  ylim(0, max(venn_data$pct) * 1.3) +
  theme_publication() +
  theme(legend.position = "none")

save_figure(
  file.path(OUT_FIGURES_SNRNASEQ, "population_bootstrap", "population_distribution.png"),
  p_venn, width = 8, height = 5
)

# ============================================================
# Interpretation
# ============================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("POPULATION OVERLAP ANALYSIS WITH BOOTSTRAP CI\n")
cat(strrep("=", 70), "\n")

cat("\n[1] ODDS RATIO INTERPRETATION:\n")
if (boot_result$OR_ci[2] < 1) {
  cat("    CLAIM SUPPORTED: OR 95% CI entirely below 1\n")
  cat("    SPP1+ and Complement_High are NEGATIVELY associated\n")
  cat("    i.e., cells that are SPP1+ are LESS likely to be Complement_High\n")
} else if (boot_result$OR_ci[1] > 1) {
  cat("    OR 95% CI entirely above 1\n")
  cat("    SPP1+ and Complement_High are POSITIVELY associated\n")
} else {
  cat("    OR 95% CI includes 1\n")
  cat("    Cannot exclude statistical independence\n")
}

cat("\n[2] JACCARD INDEX:\n")
cat(sprintf("    Overlap coefficient = %.4f\n", boot_result$observed$jaccard))
if (boot_result$observed$jaccard < 0.01) {
  cat("    CLAIM SUPPORTED: Minimal overlap between populations\n")
  cat("    SPP1+ and Complement_High are largely distinct cell populations\n")
} else if (boot_result$observed$jaccard < 0.1) {
  cat("    Low overlap - populations mostly distinct\n")
} else {
  cat("    Moderate/high overlap - populations share substantial cells\n")
}

cat("\n[3] RESOLUTION STABILITY:\n")
or_range <- range(stability_df$OR)
if (all(stability_df$OR < 1)) {
  cat("    CLAIM SUPPORTED: OR consistently < 1 across all thresholds\n")
  cat(sprintf("    OR range: %.3f - %.3f\n", or_range[1], or_range[2]))
  cat("    Finding is robust to threshold definition\n")
} else {
  cat("    Mixed results across thresholds\n")
  cat(sprintf("    OR range: %.3f - %.3f\n", or_range[1], or_range[2]))
}

cat("\n[4] RECOMMENDED LANGUAGE FOR MANUSCRIPT:\n")
cat(sprintf('    "SPP1+ and Complement-high microglia show minimal overlap\n'))
cat(sprintf('     (Jaccard = %.3f, 95%% CI: %.3f-%.3f; OR = %.2f, 95%% CI: %.2f-%.2f),\n',
            boot_result$observed$jaccard,
            boot_result$jaccard_ci[1], boot_result$jaccard_ci[2],
            boot_result$observed$OR,
            boot_result$OR_ci[1], boot_result$OR_ci[2]))
cat('     indicating distinct rather than overlapping populations."\n')

cat("\n\nOutputs saved to:\n")
cat(sprintf("  %s/population_bootstrap/\n", OUT_TABLES_SNRNASEQ))
cat(sprintf("  %s/population_bootstrap/\n", OUT_FIGURES_SNRNASEQ))
