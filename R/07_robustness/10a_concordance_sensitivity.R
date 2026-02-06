# 10a_concordance_sensitivity.R
# Concordance sensitivity analysis at relaxed thresholds
# Tests whether 50-gene set is floor or ceiling for conserved signature
# Input: output/tables/comparison/cross_platform_validation.csv
# Output: output/tables/comparison/concordance_sensitivity.csv

source("R/config.R")

cat("=== Concordance Sensitivity Analysis ===\n\n")

# Load the 557-gene matched dataset
cpv <- read.csv("output/tables/comparison/cross_platform_validation.csv",
                 stringsAsFactors = FALSE)
cat(sprintf("Loaded %d matched genes (polyimide implant-associated → silicon)\n\n", nrow(cpv)))

# Current strict threshold: FDR < 0.05 in BOTH platforms + concordant direction
strict <- cpv[cpv$same_dir & !is.na(cpv$padj) & cpv$padj < 0.05, ]
cat(sprintf("Strict (FDR<0.05 both, concordant): %d genes\n", nrow(strict)))

# Define sensitivity tiers
results <- list()

# Tier 0: Current strict (FDR < 0.05 in both, concordant direction)
results[["FDR005_both"]] <- data.frame(
  tier = "FDR<0.05 both",
  polyimide_thresh = "FDR<0.05",
  silicon_thresh = "FDR<0.05",
  require_concordant = TRUE,
  n_genes = nrow(strict),
  n_up = sum(strict$category == "Implant_Up"),
  n_down = sum(strict$category == "Implant_Down")
)

# Tier 1: FDR < 0.05 in polyimide, FDR < 0.10 in silicon
t1 <- cpv[cpv$same_dir & !is.na(cpv$padj) & cpv$padj < 0.10, ]
results[["FDR005_poly_FDR010_sil"]] <- data.frame(
  tier = "FDR<0.05 poly, FDR<0.10 sil",
  polyimide_thresh = "FDR<0.05",
  silicon_thresh = "FDR<0.10",
  require_concordant = TRUE,
  n_genes = nrow(t1),
  n_up = sum(t1$category == "Implant_Up"),
  n_down = sum(t1$category == "Implant_Down")
)

# Tier 2: FDR < 0.05 in polyimide, nominal p < 0.05 in silicon
t2 <- cpv[cpv$same_dir & !is.na(cpv$pvalue) & cpv$pvalue < 0.05, ]
results[["FDR005_poly_p005_sil"]] <- data.frame(
  tier = "FDR<0.05 poly, p<0.05 sil",
  polyimide_thresh = "FDR<0.05",
  silicon_thresh = "p<0.05",
  require_concordant = TRUE,
  n_genes = nrow(t2),
  n_up = sum(t2$category == "Implant_Up"),
  n_down = sum(t2$category == "Implant_Down")
)

# Tier 3: FDR < 0.05 in polyimide, nominal p < 0.10 in silicon
t3 <- cpv[cpv$same_dir & !is.na(cpv$pvalue) & cpv$pvalue < 0.10, ]
results[["FDR005_poly_p010_sil"]] <- data.frame(
  tier = "FDR<0.05 poly, p<0.10 sil",
  polyimide_thresh = "FDR<0.05",
  silicon_thresh = "p<0.10",
  require_concordant = TRUE,
  n_genes = nrow(t3),
  n_up = sum(t3$category == "Implant_Up"),
  n_down = sum(t3$category == "Implant_Down")
)

# Tier 4: FDR < 0.01 in both (more stringent)
t4 <- cpv[cpv$same_dir & !is.na(cpv$padj) & cpv$padj < 0.01 &
          cpv$fdr_impl_ctrl < 0.01, ]
results[["FDR001_both"]] <- data.frame(
  tier = "FDR<0.01 both",
  polyimide_thresh = "FDR<0.01",
  silicon_thresh = "FDR<0.01",
  require_concordant = TRUE,
  n_genes = nrow(t4),
  n_up = sum(t4$category == "Implant_Up"),
  n_down = sum(t4$category == "Implant_Down")
)

# Tier 5: Direction concordant only (no significance filter in silicon)
t5 <- cpv[cpv$same_dir, ]
results[["concordant_direction_only"]] <- data.frame(
  tier = "Concordant direction only",
  polyimide_thresh = "FDR<0.05",
  silicon_thresh = "none",
  require_concordant = TRUE,
  n_genes = nrow(t5),
  n_up = sum(t5$category == "Implant_Up"),
  n_down = sum(t5$category == "Implant_Down")
)

# Combine results
sensitivity_df <- do.call(rbind, results)
rownames(sensitivity_df) <- NULL

# Add percentage of matched genes
sensitivity_df$pct_of_matched <- round(100 * sensitivity_df$n_genes / nrow(cpv), 1)
sensitivity_df$pct_upregulated <- round(100 * sensitivity_df$n_up / sensitivity_df$n_genes, 1)

# Print
cat("\n=== CONCORDANCE SENSITIVITY ANALYSIS ===\n\n")
for (i in seq_len(nrow(sensitivity_df))) {
  row <- sensitivity_df[i, ]
  cat(sprintf("%-35s  %3d genes (%4.1f%% of 557)  [%d up, %d down]\n",
              row$tier, row$n_genes, row$pct_of_matched, row$n_up, row$n_down))
}

# Additional analysis: what genes are gained at relaxed thresholds?
new_at_relaxed <- setdiff(t2$gene_upper, strict$gene_upper)
if (length(new_at_relaxed) > 0) {
  cat(sprintf("\n%d genes gained at relaxed threshold (FDR poly, p<0.05 sil):\n", length(new_at_relaxed)))
  gained <- cpv[cpv$gene_upper %in% new_at_relaxed, c("gene_upper", "lfc_impl_ctrl",
                                                        "log2FoldChange", "pvalue", "padj")]
  gained <- gained[order(gained$pvalue), ]
  cat(sprintf("  Top 10 by silicon p-value:\n"))
  print(head(gained, 10))
}

# Check genome-wide direction concordance at each tier
cat("\n\n=== GENOME-WIDE CONTEXT ===\n")
cat(sprintf("Total matched genes: %d\n", nrow(cpv)))
cat(sprintf("Direction concordant (any): %d (%.1f%%)\n",
            sum(cpv$same_dir, na.rm = TRUE),
            100 * mean(cpv$same_dir, na.rm = TRUE)))
cat(sprintf("Direction concordant + silicon p<0.05: %d\n", nrow(t2)))
cat(sprintf("Direction concordant + silicon FDR<0.05: %d (current threshold)\n", nrow(strict)))

# Save
write.csv(sensitivity_df, "output/tables/comparison/concordance_sensitivity.csv",
          row.names = FALSE)

# Also save the gained genes list
if (length(new_at_relaxed) > 0) {
  gained_full <- cpv[cpv$gene_upper %in% new_at_relaxed, ]
  gained_full <- gained_full[order(gained_full$pvalue), ]
  write.csv(gained_full, "output/tables/comparison/concordance_relaxed_genes.csv",
            row.names = FALSE)
}

cat("\n\nOutputs saved:\n")
cat("  output/tables/comparison/concordance_sensitivity.csv\n")
cat("  output/tables/comparison/concordance_relaxed_genes.csv\n")
cat("\nDone.\n")
