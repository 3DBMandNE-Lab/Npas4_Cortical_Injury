# 05m_spatial_subtype_validation.R
# Spatial validation of neuronal subtype and oligodendrocyte findings
# Key questions:
# 1. Do excitatory neuron markers show stronger spatial depletion than inhibitory?
# 2. Does Qk show spatial depletion near electrodes?

library(Seurat)
library(dplyr)
library(ggplot2)
source("R/config.R")

cat("=== Spatial Validation of Subtype Findings ===\n\n")

out_dir <- file.path(OUT_TABLES_SPATIAL, "subtype_validation")
fig_dir <- file.path(OUT_FIGURES_MANUSCRIPT)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load spatial data
spatial_list <- readRDS("data/processed/spatial_seurat_list.RDS")
cat(sprintf("Loaded %d Visium samples\n", length(spatial_list)))

# Define signatures based on snRNA-seq findings
signatures <- list(
  # Excitatory - genes significantly downregulated in snRNA-seq
  Excitatory_Sig = c("Camk2a", "Grin2a", "Slc17a7", "Arc", "Snap25", "Syt1"),

  # Inhibitory - markers (these had 0 sig genes in snRNA-seq)
  Inhibitory_Sig = c("Gad1", "Gad2", "Slc32a1", "Pvalb", "Sst"),

  # Premyelinating/Qk - the significant finding
  Premyelinating_Sig = c("Qk"),  # Single gene - the significant finding

  # Mature oligodendrocyte markers
  Mature_Oligo_Sig = c("Mbp", "Plp1", "Mog", "Mag"),

  # OPC markers
  OPC_Sig = c("Pdgfra", "Cspg4"),

  # Inflammatory reference (for correlation)
  Implant_Up = c("C1qa", "C1qb", "C1qc", "C3", "Tyrobp", "Trem2", "Aif1", "Spp1")
)

# Calculate module scores
results_list <- list()

for (sample_name in names(spatial_list)) {
  cat(sprintf("\nProcessing %s...\n", sample_name))
  seu <- spatial_list[[sample_name]]

  # Check gene availability
  available_genes <- rownames(seu)

  for (sig_name in names(signatures)) {
    genes <- signatures[[sig_name]]
    genes_found <- genes[genes %in% available_genes]
    cat(sprintf("  %s: %d/%d genes found\n", sig_name, length(genes_found), length(genes)))

    if (length(genes_found) >= 1) {
      seu <- AddModuleScore(seu, features = list(genes_found),
                           name = sig_name, nbin = 24, ctrl = 100)
    }
  }

  # Extract scores
  meta <- seu@meta.data

  # Get condition
  condition <- ifelse(grepl("5B", sample_name), "Craniotomy Control",
               ifelse(grepl("4", sample_name), "Acute Stimulation",
               ifelse(grepl("7C", sample_name), "Chronic No Stim", "Chronic Stimulation")))

  # Store results
  score_cols <- grep("_Sig1$|Implant_Up1$", colnames(meta), value = TRUE)

  for (col in score_cols) {
    clean_name <- gsub("1$", "", col)
    results_list[[length(results_list) + 1]] <- data.frame(
      sample = sample_name,
      condition = condition,
      signature = clean_name,
      mean_score = mean(meta[[col]], na.rm = TRUE),
      sd_score = sd(meta[[col]], na.rm = TRUE),
      n_spots = nrow(meta)
    )
  }

  # Calculate correlations between signatures
  if ("Implant_Up1" %in% colnames(meta)) {
    for (sig_col in score_cols[score_cols != "Implant_Up1"]) {
      clean_name <- gsub("1$", "", sig_col)
      cor_test <- cor.test(meta$Implant_Up1, meta[[sig_col]], method = "pearson")

      results_list[[length(results_list) + 1]] <- data.frame(
        sample = sample_name,
        condition = condition,
        signature = paste0(clean_name, "_vs_Inflammation"),
        mean_score = cor_test$estimate,
        sd_score = NA,
        n_spots = nrow(meta),
        p_value = cor_test$p.value
      )
    }
  }

  spatial_list[[sample_name]] <- seu
}

# Combine results
results_df <- bind_rows(results_list)

# Separate scores and correlations
scores_df <- results_df %>% filter(!grepl("_vs_Inflammation", signature))
corr_df <- results_df %>% filter(grepl("_vs_Inflammation", signature))

cat("\n\n=== SPATIAL SIGNATURE SCORES ===\n")
print(scores_df %>% group_by(signature, condition) %>%
        summarize(mean = mean(mean_score), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = condition, values_from = mean))

cat("\n\n=== CORRELATIONS WITH INFLAMMATION ===\n")
corr_summary <- corr_df %>%
  mutate(signature = gsub("_vs_Inflammation", "", signature)) %>%
  group_by(signature) %>%
  summarize(
    mean_r = mean(mean_score),
    min_r = min(mean_score),
    max_r = max(mean_score),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(mean_r)

print(corr_summary)

# Key validation test: Excitatory vs Inhibitory correlation with inflammation
cat("\n\n=== KEY VALIDATION: EXCITATORY VS INHIBITORY ===\n")

exc_corr <- corr_df %>% filter(signature == "Excitatory_Sig_vs_Inflammation")
inh_corr <- corr_df %>% filter(signature == "Inhibitory_Sig_vs_Inflammation")

cat("\nExcitatory ~ Inflammation:\n")
cat(sprintf("  Mean r = %.3f (range: %.3f to %.3f)\n",
            mean(exc_corr$mean_score), min(exc_corr$mean_score), max(exc_corr$mean_score)))

cat("\nInhibitory ~ Inflammation:\n")
cat(sprintf("  Mean r = %.3f (range: %.3f to %.3f)\n",
            mean(inh_corr$mean_score), min(inh_corr$mean_score), max(inh_corr$mean_score)))

# Statistical test: are excitatory markers more negatively correlated?
if (nrow(exc_corr) >= 3 && nrow(inh_corr) >= 3) {
  test <- wilcox.test(exc_corr$mean_score, inh_corr$mean_score, paired = TRUE)
  cat(sprintf("\nPaired Wilcoxon test (Excitatory vs Inhibitory): p = %.4f\n", test$p.value))
}

# Qk validation
cat("\n\n=== QK SPATIAL VALIDATION ===\n")
qk_corr <- corr_df %>% filter(signature == "Premyelinating_Sig_vs_Inflammation")
if (nrow(qk_corr) > 0) {
  cat(sprintf("Qk ~ Inflammation: mean r = %.3f (range: %.3f to %.3f)\n",
              mean(qk_corr$mean_score), min(qk_corr$mean_score), max(qk_corr$mean_score)))
  cat("Negative correlation = Qk depleted in inflammatory zones\n")
}

# Save results
write.csv(scores_df, file.path(out_dir, "spatial_subtype_scores.csv"), row.names = FALSE)
write.csv(corr_df, file.path(out_dir, "spatial_subtype_correlations.csv"), row.names = FALSE)
write.csv(corr_summary, file.path(out_dir, "spatial_correlation_summary.csv"), row.names = FALSE)

cat(sprintf("\nResults saved to: %s/\n", out_dir))

# Create validation figure
cat("\n\nGenerating validation figure...\n")

corr_plot_df <- corr_df %>%
  mutate(signature = gsub("_vs_Inflammation", "", signature),
         signature = gsub("_Sig", "", signature))

p <- ggplot(corr_plot_df, aes(x = reorder(signature, mean_score), y = mean_score)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  geom_point(aes(color = condition), size = 2) +
  coord_flip() +
  labs(
    title = "Spatial Correlation with Inflammatory Signature",
    subtitle = "Negative = depleted in inflammatory zones",
    x = "Signature",
    y = "Pearson r with Implant_Up"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "spatial_subtype_validation.png"), p,
       width = 8, height = 6, dpi = 300)

cat("\nFigure saved: spatial_subtype_validation.png\n")
cat("\nDone.\n")
