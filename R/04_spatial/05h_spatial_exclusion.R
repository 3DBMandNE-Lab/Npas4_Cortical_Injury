# 05h_spatial_exclusion.R
# Spatial exclusion/overlap analysis between enriched domains
# Input: data/processed/spatial_scored.RDS
# Output: output/tables/spatial/domain_overlap_matrix.csv, domain_exclusion_stats.csv

library(dplyr)
library(tidyr)
source("R/config.R")

seurat_list <- readRDS(file.path(DATA_PROCESSED, "spatial_scored.RDS"))

signatures <- c("Implant_Up", "DAM", "Neuronal", "Complement")

# Function to calculate domain overlap
calc_domain_overlap <- function(meta, thresh = 0.75) {
  avail <- intersect(signatures, colnames(meta))
  if (length(avail) < 2) return(NULL)

  n_spots <- nrow(meta)

  # Define "enriched" domains as top 25% for each signature
  # For Neuronal, "enriched" means low activity (bottom 25%)
  domains <- list()
  for (sig in avail) {
    if (sig == "Neuronal") {
      # Low neuronal = "enriched" in dysfunction
      domains[[sig]] <- meta[[sig]] < quantile(meta[[sig]], 1 - thresh, na.rm = TRUE)
    } else {
      # High inflammation = enriched
      domains[[sig]] <- meta[[sig]] > quantile(meta[[sig]], thresh, na.rm = TRUE)
    }
  }

  # Calculate pairwise overlap statistics
  results <- list()
  for (i in 1:length(avail)) {
    for (j in 1:length(avail)) {
      sig1 <- avail[i]
      sig2 <- avail[j]

      d1 <- domains[[sig1]]
      d2 <- domains[[sig2]]

      # Counts
      n_sig1 <- sum(d1, na.rm = TRUE)
      n_sig2 <- sum(d2, na.rm = TRUE)
      n_both <- sum(d1 & d2, na.rm = TRUE)
      n_sig1_only <- sum(d1 & !d2, na.rm = TRUE)
      n_sig2_only <- sum(!d1 & d2, na.rm = TRUE)
      n_neither <- sum(!d1 & !d2, na.rm = TRUE)

      # Jaccard index: intersection / union
      union_size <- n_sig1 + n_sig2 - n_both
      jaccard <- if (union_size > 0) n_both / union_size else NA

      # Overlap coefficient: intersection / min(size1, size2)
      min_size <- min(n_sig1, n_sig2)
      overlap_coef <- if (min_size > 0) n_both / min_size else NA

      # Expected overlap under independence
      expected_both <- (n_sig1 / n_spots) * (n_sig2 / n_spots) * n_spots
      fold_enrichment <- if (expected_both > 0) n_both / expected_both else NA

      # Fisher's exact test for enrichment/depletion
      contingency <- matrix(c(n_both, n_sig1_only, n_sig2_only, n_neither), nrow = 2)
      fisher_p <- tryCatch(fisher.test(contingency)$p.value, error = function(e) NA)

      results[[paste(sig1, sig2, sep = "_")]] <- data.frame(
        domain1 = sig1,
        domain2 = sig2,
        n_domain1 = n_sig1,
        n_domain2 = n_sig2,
        n_overlap = n_both,
        n_domain1_only = n_sig1_only,
        n_domain2_only = n_sig2_only,
        pct_domain1_in_overlap = 100 * n_both / n_sig1,
        pct_domain2_in_overlap = 100 * n_both / n_sig2,
        jaccard_index = jaccard,
        overlap_coefficient = overlap_coef,
        expected_overlap = expected_both,
        fold_enrichment = fold_enrichment,
        fisher_pvalue = fisher_p
      )
    }
  }

  bind_rows(results)
}

# Process all samples
all_results <- list()

for (sid in names(seurat_list)) {
  cat(sprintf("Processing: %s\n", sid))
  meta <- seurat_list[[sid]]@meta.data

  overlap_df <- calc_domain_overlap(meta)
  if (is.null(overlap_df)) next

  overlap_df$sample <- sid
  overlap_df$condition <- meta$condition[1]
  all_results[[sid]] <- overlap_df
}

results_df <- bind_rows(all_results)

# Create summary matrix for each sample
cat("\n=== Domain Overlap Summary ===\n")

for (sid in names(seurat_list)) {
  sub <- results_df[results_df$sample == sid, ]
  if (nrow(sub) == 0) next

  cat(sprintf("\n%s (%s):\n", sid, sub$condition[1]))

  # Print Jaccard matrix
  jaccard_mat <- sub %>%
    select(domain1, domain2, jaccard_index) %>%
    pivot_wider(names_from = domain2, values_from = jaccard_index)
  print(as.data.frame(jaccard_mat))
}

# Summary of exclusion patterns
cat("\n=== Exclusion Patterns (Fold Enrichment < 0.5) ===\n")
exclusions <- results_df %>%
  filter(domain1 != domain2, fold_enrichment < 0.5) %>%
  select(sample, condition, domain1, domain2, fold_enrichment, fisher_pvalue)
print(exclusions)

cat("\n=== Strong Colocalization (Fold Enrichment > 2) ===\n")
colocalizations <- results_df %>%
  filter(domain1 != domain2, fold_enrichment > 2) %>%
  select(sample, condition, domain1, domain2, fold_enrichment, fisher_pvalue)
print(colocalizations)

# Save results
write.csv(results_df, file.path(OUT_TABLES_SPATIAL, "domain_overlap_matrix.csv"), row.names = FALSE)

# Summary statistics across samples
summary_stats <- results_df %>%
  filter(domain1 != domain2) %>%
  group_by(domain1, domain2) %>%
  summarize(
    mean_jaccard = mean(jaccard_index, na.rm = TRUE),
    sd_jaccard = sd(jaccard_index, na.rm = TRUE),
    mean_fold_enrichment = mean(fold_enrichment, na.rm = TRUE),
    sd_fold_enrichment = sd(fold_enrichment, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

print(summary_stats)
write.csv(summary_stats, file.path(OUT_TABLES_SPATIAL, "domain_exclusion_summary.csv"), row.names = FALSE)

cat(sprintf("\nSaved: %s/domain_overlap_matrix.csv\n", OUT_TABLES_SPATIAL))
cat(sprintf("Saved: %s/domain_exclusion_summary.csv\n", OUT_TABLES_SPATIAL))
