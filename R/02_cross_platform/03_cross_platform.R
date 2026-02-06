# 03_cross_platform.R
# Validate polyimide signature in silicon data
# Input: comparison/implant_specific.csv, deg/silicon_deseq2_results.csv
# Output: comparison/cross_platform_validation.csv, validated_genes.csv

source("R/config.R")

# Load results
implant <- read.csv(file.path(OUT_TABLES_COMPARISON, "implant_specific.csv"))
silicon <- read.csv(file.path(OUT_TABLES_DEG, "silicon_deseq2_results.csv"))

cat(sprintf("Implant-specific genes: %d\n", nrow(implant)))
cat(sprintf("Silicon genes: %d\n", nrow(silicon)))

# Match genes (case-insensitive)
implant$gene_upper <- toupper(implant$gene)
silicon$gene_upper <- toupper(silicon$gene)

matched <- merge(implant, silicon, by = "gene_upper", suffixes = c("_poly", "_sil"))
cat(sprintf("Matched genes: %d (%.1f%%)\n", nrow(matched), 100 * nrow(matched) / nrow(implant)))

# Direction concordance
matched$same_dir <- sign(matched$lfc_impl_ctrl) == sign(matched$log2FoldChange)
n_concordant <- sum(matched$same_dir, na.rm = TRUE)
n_total <- sum(!is.na(matched$same_dir))
pct_concordant <- 100 * n_concordant / n_total

cat(sprintf("\nDirection concordance: %d/%d (%.1f%%)\n", n_concordant, n_total, pct_concordant))

# Binomial test
binom_p <- binom.test(n_concordant, n_total, p = 0.5, alternative = "greater")$p.value
cat(sprintf("Binomial test p-value: %.2e\n", binom_p))

# Effect size correlation
cor_test <- cor.test(matched$lfc_impl_ctrl, matched$log2FoldChange, method = "spearman")
cat(sprintf("Spearman rho: %.3f (p = %.2e)\n", cor_test$estimate, cor_test$p.value))

# Validated genes (significant in both, same direction)
validated <- matched[matched$same_dir & matched$significant, ]
cat(sprintf("\nValidated genes (sig in silicon + concordant): %d\n", nrow(validated)))

# Save
write.csv(matched, file.path(OUT_TABLES_COMPARISON, "cross_platform_validation.csv"), row.names = FALSE)
write.csv(validated, file.path(OUT_TABLES_COMPARISON, "validated_genes.csv"), row.names = FALSE)

cat("\nSaved to: output/tables/comparison/\n")
