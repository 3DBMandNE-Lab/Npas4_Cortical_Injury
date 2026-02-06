# 10b_bootstrap_ci_check.R
# Bootstrap CI method verification for Table 1
# Compares: (a) current delta method, (b) percentile bootstrap, (c) BCa bootstrap
# Input: output/tables/snrnaseq/signature_scores_summary.csv
# Output: output/tables/snrnaseq/bootstrap_ci_comparison.csv

# Compares delta method, percentile bootstrap, and BCa bootstrap CI methods

source("R/config.R")

cat("=== Bootstrap CI Method Check (Table 1) ===\n\n")

# Load precomputed scores
scores <- read.csv("output/tables/snrnaseq/signature_scores_summary.csv",
                    stringsAsFactors = FALSE)

# Five key programs from Table 1
key_programs <- c("Calcium_Signaling", "Survival", "Synaptic_Core",
                  "Glutamatergic", "Activity_IEGs")

# Filter to Implant and Control
impl <- scores[scores$Condition == "Implant" & scores$signature %in% key_programs, ]
ctrl <- scores[scores$Condition == "Control" & scores$signature %in% key_programs, ]

cat("Method A: Delta method CIs (current implementation in fig5b code)\n")
cat("Method B: Percentile bootstrap (as claimed in manuscript Methods)\n")
cat("Method C: BCa bootstrap (bias-corrected and accelerated)\n\n")

# ============================================================
# Method A: Delta method (current)
# ============================================================
delta_results <- merge(impl, ctrl, by = "signature", suffixes = c("_impl", "_ctrl"))
delta_results$pct_control <- 100 * delta_results$mean_score_impl / delta_results$mean_score_ctrl
delta_results$se_impl <- delta_results$sd_score_impl / sqrt(delta_results$n_cells_impl)
delta_results$se_ctrl <- delta_results$sd_score_ctrl / sqrt(delta_results$n_cells_ctrl)
delta_results$se_ratio <- 100 * sqrt(
  (delta_results$se_impl / delta_results$mean_score_ctrl)^2 +
  (delta_results$mean_score_impl * delta_results$se_ctrl / delta_results$mean_score_ctrl^2)^2
)
delta_results$ci_low_delta <- delta_results$pct_control - 1.96 * delta_results$se_ratio
delta_results$ci_high_delta <- delta_results$pct_control + 1.96 * delta_results$se_ratio

# ============================================================
# Methods B & C: Bootstrap on sample-level means
# The manuscript says CIs are "stratified by sample" — this means
# we should bootstrap at the sample level, not cell level
# Since we only have summary stats (not raw cell data), we simulate
# sample-level bootstrap using the available mean/SD/n statistics.
#
# However, without the raw per-sample data, we can only do a
# parametric bootstrap (assuming normality of sample means by CLT).
# ============================================================

# For a proper bootstrap we need per-sample means, not cell-level summary.
# Let's check if per-sample data exists:
sample_file <- "output/tables/snrnaseq/signature_scores_by_sample.csv"
if (!file.exists(sample_file)) {
  cat("Per-sample score file not found. Checking alternative...\n")
  # Try other possible locations
  alt_files <- list.files("output/tables/snrnaseq/", pattern = "sample",
                          full.names = TRUE, recursive = TRUE)
  if (length(alt_files) > 0) {
    cat("Found alternatives:\n")
    cat(paste("  ", alt_files, collapse = "\n"), "\n")
  }
}

# ============================================================
# Parametric bootstrap (CLT-based) as a check
# ============================================================
set.seed(42)
N_BOOT <- 10000

cat("Running parametric bootstrap (N = 10,000)...\n\n")

boot_results <- list()

for (prog in key_programs) {
  impl_row <- delta_results[delta_results$signature == prog, ]

  # Sample means from normal distribution using CLT
  # mu = mean_score, se = sd/sqrt(n)
  impl_means <- rnorm(N_BOOT, mean = impl_row$mean_score_impl,
                       sd = impl_row$sd_score_impl / sqrt(impl_row$n_cells_impl))
  ctrl_means <- rnorm(N_BOOT, mean = impl_row$mean_score_ctrl,
                       sd = impl_row$sd_score_ctrl / sqrt(impl_row$n_cells_ctrl))

  # Ratio: % of control
  boot_ratios <- 100 * impl_means / ctrl_means

  # Remove infinite/NaN (ctrl_mean near zero)
  boot_ratios <- boot_ratios[is.finite(boot_ratios)]

  # Percentile CI
  pct_ci <- quantile(boot_ratios, c(0.025, 0.975))

  # BCa requires jackknife — with parametric bootstrap, we approximate:
  # z0 = bias correction
  observed_ratio <- impl_row$pct_control
  z0 <- qnorm(mean(boot_ratios < observed_ratio))

  # Acceleration (jackknife approximation)
  # For parametric bootstrap, acceleration is approximately 0
  # (the parametric model has no skewness in the mean estimator by CLT)
  # So BCa ≈ bias-corrected percentile with z0 adjustment
  a_hat <- 0  # no acceleration for parametric
  alpha <- c(0.025, 0.975)
  z_alpha <- qnorm(alpha)
  adjusted_alpha <- pnorm(z0 + (z0 + z_alpha) / (1 - a_hat * (z0 + z_alpha)))
  bca_ci <- quantile(boot_ratios, adjusted_alpha)

  boot_results[[prog]] <- data.frame(
    signature = prog,
    pct_control = round(observed_ratio, 1),
    # Delta method
    ci_low_delta = round(impl_row$ci_low_delta, 1),
    ci_high_delta = round(impl_row$ci_high_delta, 1),
    # Percentile bootstrap
    ci_low_pct = round(pct_ci[1], 1),
    ci_high_pct = round(pct_ci[2], 1),
    # BCa bootstrap
    ci_low_bca = round(bca_ci[1], 1),
    ci_high_bca = round(bca_ci[2], 1),
    # Diagnostics
    boot_mean = round(mean(boot_ratios), 1),
    boot_median = round(median(boot_ratios), 1),
    z0_bias = round(z0, 4),
    n_valid_boot = length(boot_ratios),
    stringsAsFactors = FALSE
  )
}

comparison <- do.call(rbind, boot_results)
rownames(comparison) <- NULL

# ============================================================
# Compare with manuscript Table 1 values
# ============================================================
table1 <- data.frame(
  signature = c("Calcium_Signaling", "Survival", "Synaptic_Core",
                "Glutamatergic", "Activity_IEGs"),
  manuscript_pct = c(22, 15, 67, 58, 13),
  manuscript_ci_low = c(19, 12, 64, 55, 9),
  manuscript_ci_high = c(25, 18, 70, 61, 17)
)

comparison <- merge(comparison, table1, by = "signature")

cat("\n=== COMPARISON: Table 1 CIs Across Methods ===\n\n")
cat(sprintf("%-20s  %s  |  %-15s  |  %-15s  |  %-15s  |  %-15s\n",
            "Program", "% Ctrl", "Manuscript", "Delta", "Percentile", "BCa"))
cat(paste(rep("-", 110), collapse = ""), "\n")

for (i in seq_len(nrow(comparison))) {
  r <- comparison[i, ]
  cat(sprintf("%-20s  %5.1f  |  %4.0f (%3.0f-%3.0f)  |  %5.1f (%5.1f-%5.1f)  |  %5.1f (%5.1f-%5.1f)  |  %5.1f (%5.1f-%5.1f)\n",
              r$signature, r$pct_control,
              r$manuscript_pct, r$manuscript_ci_low, r$manuscript_ci_high,
              r$pct_control, r$ci_low_delta, r$ci_high_delta,
              r$pct_control, r$ci_low_pct, r$ci_high_pct,
              r$pct_control, r$ci_low_bca, r$ci_high_bca))
}

# ============================================================
# Assessment
# ============================================================
cat("\n\n=== ASSESSMENT ===\n\n")

cat("1. MANUSCRIPT CLAIMS: '10,000 iterations, percentile method'\n")
cat("   ACTUAL CODE: Delta method (normal approximation) in fig5b\n")
cat("   DISCREPANCY: Yes — Methods do not match implementation\n\n")

cat("2. PRACTICAL IMPACT:\n")
# Check if CIs differ substantially
max_diff <- max(abs(comparison$ci_low_delta - comparison$ci_low_pct),
                abs(comparison$ci_high_delta - comparison$ci_high_pct))
cat(sprintf("   Max difference between delta and percentile: %.1f percentage points\n", max_diff))

if (max_diff < 2) {
  cat("   IMPACT: Minimal — all methods agree within ~2 pp\n")
  cat("   RECOMMENDATION: Update Methods to state 'delta method (normal approximation)'\n")
  cat("   or actually implement the bootstrap. Both approaches are defensible.\n")
} else {
  cat("   IMPACT: Moderate — methods diverge meaningfully\n")
  cat("   RECOMMENDATION: Implement actual sample-level bootstrap from raw data\n")
}

cat("\n3. BCa vs PERCENTILE:\n")
max_bca_diff <- max(abs(comparison$ci_low_pct - comparison$ci_low_bca),
                    abs(comparison$ci_high_pct - comparison$ci_high_bca))
cat(sprintf("   Max difference: %.1f percentage points\n", max_bca_diff))
cat("   For parametric bootstrap, BCa ≈ percentile (acceleration = 0 by CLT)\n")
cat("   BCa would differ more with non-parametric sample-level bootstrap\n")

cat("\n4. NOTE ON Activity_IEGs:\n")
cat("   Both Implant and Control means are near zero (~-0.01 and ~-0.0005)\n")
cat("   The 'ratio' is unstable: -0.01 / -0.0005 = 18.2 (not 13%)\n")
cat("   The manuscript value (13%) appears to be computed differently\n")
cat("   This program requires careful handling due to near-zero denominator\n")

# Save
write.csv(comparison, "output/tables/snrnaseq/bootstrap_ci_comparison.csv",
          row.names = FALSE)

cat("\n\nOutput saved: output/tables/snrnaseq/bootstrap_ci_comparison.csv\n")
cat("Done.\n")
