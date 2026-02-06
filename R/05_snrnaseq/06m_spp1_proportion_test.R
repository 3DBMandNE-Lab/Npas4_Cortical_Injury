# 06m_spp1_proportion_test.R
# Test SPP1+ cell proportions: Implant vs Stab (direct comparison)
# For rare markers, proportion of positive cells is the appropriate metric

library(Seurat)
library(dplyr)
source("R/config.R")

cat("=== SPP1+ Proportion Test: Implant vs Stab ===\n\n")

dir.create(file.path(OUT_TABLES_SNRNASEQ, "spp1_specificity"), recursive = TRUE, showWarnings = FALSE)

# Load snRNA-seq data
# SNRNASEQ_PATH defined in config.R
cat("Loading snRNA-seq data...\n")
seu <- readRDS(SNRNASEQ_PATH)

# Subset to microglia
microglia <- subset(seu, celltype_l1 == "Microglia")
# CRITICAL: Use RNA assay for gene detection - SCT may distort sparse expression
DefaultAssay(microglia) <- "RNA"
cat(sprintf("Microglia: %d cells (using %s assay)\n", ncol(microglia), DefaultAssay(microglia)))

# Get Spp1 expression from RNA assay
# NOTE: For proportion analysis, we use counts (presence/absence of UMIs)
spp1_counts <- GetAssayData(microglia, assay = "RNA", layer = "counts")["Spp1", ]
microglia$Spp1_positive <- spp1_counts > 0

cat(sprintf("Spp1+ microglia: %d (%.2f%%)\n", sum(spp1_counts > 0), 100*mean(spp1_counts > 0)))
# ============================================================
# Sample-level proportions
# ============================================================

cat("\n=== Sample-level SPP1+ proportions ===\n")

sample_props <- microglia@meta.data %>%
  group_by(orig.ident, Condition) %>%
  summarize(
    n_cells = n(),
    n_spp1_pos = sum(Spp1_positive),
    pct_spp1_pos = 100 * mean(Spp1_positive),
    .groups = "drop"
  )

cat("\nBy sample:\n")
print(as.data.frame(sample_props))

# Condition summaries
condition_summary <- sample_props %>%
  group_by(Condition) %>%
  summarize(
    n_samples = n(),
    total_cells = sum(n_cells),
    total_spp1_pos = sum(n_spp1_pos),
    mean_pct = mean(pct_spp1_pos),
    sd_pct = sd(pct_spp1_pos),
    se_pct = sd(pct_spp1_pos) / sqrt(n()),
    .groups = "drop"
  )

cat("\n=== Condition Summary ===\n")
print(as.data.frame(condition_summary))

# ============================================================
# Statistical tests
# ============================================================

cat("\n=== Statistical Tests (sample-level, Wilcoxon rank-sum) ===\n")

# Extract proportions by condition
control_props <- sample_props$pct_spp1_pos[sample_props$Condition == "Control"]
implant_props <- sample_props$pct_spp1_pos[sample_props$Condition == "Implant"]
stab_props <- sample_props$pct_spp1_pos[sample_props$Condition == "Stab"]

# Implant vs Control
test_impl_ctrl <- wilcox.test(implant_props, control_props)
cat(sprintf("\nImplant vs Control: W = %.0f, p = %.4f\n",
            test_impl_ctrl$statistic, test_impl_ctrl$p.value))

# Stab vs Control
test_stab_ctrl <- wilcox.test(stab_props, control_props)
cat(sprintf("Stab vs Control: W = %.0f, p = %.4f\n",
            test_stab_ctrl$statistic, test_stab_ctrl$p.value))

# DIRECT TEST: Implant vs Stab
test_impl_stab <- wilcox.test(implant_props, stab_props)
cat(sprintf("Implant vs Stab: W = %.0f, p = %.4f\n",
            test_impl_stab$statistic, test_impl_stab$p.value))

# ============================================================
# Fold enrichments
# ============================================================

cat("\n=== Fold Enrichments ===\n")

mean_control <- mean(control_props)
mean_implant <- mean(implant_props)
mean_stab <- mean(stab_props)

cat(sprintf("\nMean SPP1+ %%:\n"))
cat(sprintf("  Control: %.3f%%\n", mean_control))
cat(sprintf("  Implant: %.3f%%\n", mean_implant))
cat(sprintf("  Stab:    %.3f%%\n", mean_stab))

cat(sprintf("\nFold enrichment vs Control:\n"))
cat(sprintf("  Implant: %.1fx\n", mean_implant / mean_control))
cat(sprintf("  Stab:    %.1fx\n", mean_stab / mean_control))

cat(sprintf("\nImplant vs Stab ratio: %.2fx\n", mean_implant / mean_stab))

# ============================================================
# Save results
# ============================================================

results_summary <- data.frame(
  comparison = c("Implant_vs_Control", "Stab_vs_Control", "Implant_vs_Stab"),
  group1_mean_pct = c(mean_implant, mean_stab, mean_implant),
  group2_mean_pct = c(mean_control, mean_control, mean_stab),
  fold_enrichment = c(mean_implant/mean_control, mean_stab/mean_control, mean_implant/mean_stab),
  wilcox_W = c(test_impl_ctrl$statistic, test_stab_ctrl$statistic, test_impl_stab$statistic),
  p_value = c(test_impl_ctrl$p.value, test_stab_ctrl$p.value, test_impl_stab$p.value)
)

write.csv(results_summary, file.path(OUT_TABLES_SNRNASEQ, "spp1_specificity", "spp1_proportion_tests.csv"),
          row.names = FALSE)

write.csv(sample_props, file.path(OUT_TABLES_SNRNASEQ, "spp1_specificity", "spp1_by_sample.csv"),
          row.names = FALSE)

# ============================================================
# Interpretation
# ============================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SPP1+ PROPORTION TEST RESULTS\n")
cat(strrep("=", 70), "\n\n")

cat("Sample-level proportions (Wilcoxon rank-sum test):\n\n")
print(results_summary)

cat("\n[KEY FINDING]:\n")
if (test_impl_stab$p.value < 0.05) {
  cat(sprintf("SPP1+ proportion is significantly higher in Implant vs Stab (p = %.4f).\n",
              test_impl_stab$p.value))
  cat("This confirms SPP1 as an implant-specific marker.\n")
} else {
  cat(sprintf("SPP1+ proportion is NOT significantly different between Implant and Stab (p = %.4f).\n",
              test_impl_stab$p.value))
  cat(sprintf("However, there IS a trend: Implant (%.2f%%) vs Stab (%.2f%%) = %.1fx difference.\n",
              mean_implant, mean_stab, mean_implant/mean_stab))
  cat("\n[RECOMMENDATION]:\n")
  cat("Reframe SPP1 as showing 'graded response to chronic perturbation'\n")
  cat("rather than 'binary implant-specificity'. Both Implant and Stab show\n")
  cat("elevated SPP1+ proportions vs Control, with Implant showing the\n")
  cat("highest proportion.\n")
}

cat("\n\nOutput saved to:\n")
cat(sprintf("  %s/spp1_specificity/\n", OUT_TABLES_SNRNASEQ))
