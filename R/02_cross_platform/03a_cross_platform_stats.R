# 03a_cross_platform_stats.R
# Full cross-platform concordance statistics
# Input: deg/*.csv, comparison/implant_specific.csv
# Output: comparison/concordance_statistics.csv, implant_signature_validation_stats.csv
#
# NOTE ON CROSS-PLATFORM COMPARISONS:
# - Absolute log2FC values are NOT directly comparable (different platforms, dynamic ranges)
# - Valid metrics: (1) Direction concordance (sign match), (2) Rank correlation (Spearman)
# - Primary metric: Direction concordance (does gene go same direction in both?)
# - Secondary metric: Rank correlation (do strongly changed genes rank similarly?)

library(dplyr)
source("R/config.R")

cat("=== Cross-Platform Concordance Analysis ===\n")
cat("NOTE: Comparing DIRECTION and RANK, not absolute effect sizes\n\n")

# Load results
silicon <- read.csv(file.path(OUT_TABLES_DEG, "silicon_deseq2_results.csv"))
polyimide <- read.csv(file.path(OUT_TABLES_DEG, "polyimide_limma_results.csv"))
implant_spec <- read.csv(file.path(OUT_TABLES_COMPARISON, "implant_specific.csv"))

cat(sprintf("Silicon genes: %d\n", nrow(silicon)))
cat(sprintf("Polyimide genes: %d\n", nrow(polyimide)))
cat(sprintf("Implant-specific genes: %d\n", nrow(implant_spec)))

# Standardize gene names
silicon$gene_upper <- toupper(silicon$gene)
polyimide$gene_upper <- toupper(polyimide$gene)
implant_spec$gene_upper <- toupper(implant_spec$gene)

# 1. ALL GENES concordance
all_matched <- merge(silicon, polyimide, by = "gene_upper", suffixes = c("_sil", "_poly"))
cat(sprintf("\nAll genes matched: %d\n", nrow(all_matched)))

# Direction concordance (PRIMARY METRIC - valid across platforms)
all_matched$same_dir <- sign(all_matched$log2FoldChange) == sign(all_matched$logFC)
dir_conc_all <- mean(all_matched$same_dir, na.rm = TRUE)
cat(sprintf("All genes - Direction concordance: %.1f%%\n", 100 * dir_conc_all))

# Rank correlation (SECONDARY METRIC - compares relative ordering, not absolute values)
# Spearman rho is based on ranks, so it's valid across platforms
rank_cor_all <- cor.test(rank(all_matched$log2FoldChange), rank(all_matched$logFC), method = "spearman")
cat(sprintf("All genes - Rank correlation (Spearman): %.3f (p = %.2e)\n",
            rank_cor_all$estimate, rank_cor_all$p.value))

# 2. IMPLANT-SPECIFIC signature validation (polyimide -> silicon)
impl_matched <- merge(implant_spec, silicon, by = "gene_upper", suffixes = c("_impl", "_sil"))
cat(sprintf("\nImplant-specific matched in silicon: %d / %d (%.1f%%)\n",
            nrow(impl_matched), nrow(implant_spec), 100 * nrow(impl_matched) / nrow(implant_spec)))

# Direction concordance for implant-specific (PRIMARY METRIC)
impl_matched$same_dir <- sign(impl_matched$lfc_impl_ctrl) == sign(impl_matched$log2FoldChange)
n_concordant <- sum(impl_matched$same_dir, na.rm = TRUE)
n_total <- sum(!is.na(impl_matched$same_dir))
dir_conc_impl <- n_concordant / n_total

cat(sprintf("Implant-specific - Direction concordance: %d/%d (%.1f%%)\n",
            n_concordant, n_total, 100 * dir_conc_impl))

# Binomial test: Is concordance better than chance (50%)?
binom_p <- binom.test(n_concordant, n_total, p = 0.5, alternative = "greater")$p.value
cat(sprintf("Binomial test (H0: concordance = 50%%): p = %.2e\n", binom_p))

# Rank correlation for implant-specific (SECONDARY METRIC)
rank_cor_impl <- cor.test(rank(impl_matched$lfc_impl_ctrl), rank(impl_matched$log2FoldChange),
                          method = "spearman")
cat(sprintf("Implant-specific - Rank correlation (Spearman): %.3f (p = %.2e)\n",
            rank_cor_impl$estimate, rank_cor_impl$p.value))

# Validated genes (significant in both platforms, same direction)
impl_matched$sig_sil <- impl_matched$padj < FDR_THRESH
validated <- impl_matched[impl_matched$same_dir & impl_matched$sig_sil & !is.na(impl_matched$sig_sil), ]
cat(sprintf("\nValidated genes (FDR < 0.05 in silicon + concordant direction): %d\n", nrow(validated)))

# Save concordance statistics
conc_stats <- data.frame(
  metric = c("all_genes_direction_concordance",
             "all_genes_rank_correlation",
             "implant_specific_matched",
             "implant_specific_direction_concordance",
             "implant_specific_rank_correlation",
             "implant_specific_validated"),
  value = c(dir_conc_all,
            rank_cor_all$estimate,
            nrow(impl_matched),
            dir_conc_impl,
            rank_cor_impl$estimate,
            nrow(validated)),
  p_value = c(NA,
              rank_cor_all$p.value,
              NA,
              binom_p,
              rank_cor_impl$p.value,
              NA),
  note = c("Fraction with same sign",
           "Spearman rho on ranks (platform-independent)",
           "Genes found in both datasets",
           "Fraction with same sign (binomial p-value)",
           "Spearman rho on ranks (platform-independent)",
           "FDR < 0.05 in both + same direction")
)

write.csv(conc_stats, file.path(OUT_TABLES_COMPARISON, "concordance_statistics.csv"), row.names = FALSE)

# Save implant signature validation details
impl_validation <- data.frame(
  n_implant_specific = nrow(implant_spec),
  n_matched_silicon = nrow(impl_matched),
  pct_matched = 100 * nrow(impl_matched) / nrow(implant_spec),
  n_direction_concordant = n_concordant,
  pct_concordant = 100 * dir_conc_impl,
  binomial_pvalue = binom_p,
  rank_correlation = rank_cor_impl$estimate,
  rank_correlation_pvalue = rank_cor_impl$p.value,
  n_validated = nrow(validated)
)

write.csv(impl_validation, file.path(OUT_TABLES_COMPARISON, "implant_signature_validation_stats.csv"), row.names = FALSE)
write.csv(validated, file.path(OUT_TABLES_COMPARISON, "validated_genes.csv"), row.names = FALSE)

# Save cross-platform validation table (for plotting)
write.csv(impl_matched, file.path(OUT_TABLES_COMPARISON, "cross_platform_validation.csv"), row.names = FALSE)

cat(sprintf("\nSaved: %s/concordance_statistics.csv\n", OUT_TABLES_COMPARISON))
cat(sprintf("Saved: %s/implant_signature_validation_stats.csv\n", OUT_TABLES_COMPARISON))
cat(sprintf("Saved: %s/validated_genes.csv (%d genes)\n", OUT_TABLES_COMPARISON, nrow(validated)))
