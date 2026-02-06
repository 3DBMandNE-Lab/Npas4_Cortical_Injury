# 06f_snrnaseq_signatures.R
# Full neuronal signature analysis with module scores
# Input: Seurat object (path set by SNRNASEQ_PATH in config.R)
# Output: tables/snrnaseq/signature_*.csv, figures/snrnaseq/signature_*.png

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/config.R")

SNRNASEQ_SOURCE <- SNRNASEQ_PATH  # defined in config.R

cat("Loading full snRNA-seq dataset...\n")
seu <- readRDS(SNRNASEQ_SOURCE)
cat(sprintf("Loaded: %d cells, %d genes\n", ncol(seu), nrow(seu)))

# Use RNA assay for module scoring (more stable than SCT)
DefaultAssay(seu) <- "RNA"
cat(sprintf("Using assay: %s\n", DefaultAssay(seu)))

# Define all gene signatures from supplementary methods
signatures <- list(
  Activity_IEGs = c("Npas4", "Arc", "Fos", "Egr1", "Nr4a1", "Nr4a3", "Junb", "Fosb"),
  Synaptic_Core = c("Rbfox3", "Gria1", "Gria2", "Shank2", "Syt1", "Homer1", "Dlg4", "Snap25"),
  Catastrophic = c("Rbfox3", "Shank2", "Nefl", "Nefm", "Nefh"),
  Chronic = c("Gria1", "Homer1", "Dlg4"),
  Preserved = c("Syt1", "Snap25", "Map2", "Tubb3"),
  Glutamatergic = c("Slc17a7", "Slc17a6", "Grin1", "Grin2a", "Grin2b", "Gria1", "Gria2"),
  GABAergic = c("Gad1", "Gad2", "Slc32a1", "Pvalb", "Cck", "Sst", "Vip"),
  Calcium_Signaling = c("Camk2a", "Camk2b", "Calm1", "Cacna1c", "Cacna1d"),
  Mitochondrial = c("Ndufa1", "Cox6c", "Atp5g1", "Sod2"),
  Apoptosis = c("Bax", "Bcl2", "Casp3", "Casp9", "Apaf1"),
  Survival = c("Bdnf", "Ntrk2", "Akt1", "Pik3ca")
)

available_genes <- rownames(seu)

# Check gene availability
cat("\nGene availability per signature:\n")
for (sig_name in names(signatures)) {
  genes <- signatures[[sig_name]]
  present <- intersect(genes, available_genes)
  cat(sprintf("  %s: %d/%d genes present\n", sig_name, length(present), length(genes)))
}

# Score each signature using AddModuleScore
cat("\nScoring signatures...\n")
scored_sigs <- c()
for (sig_name in names(signatures)) {
  genes <- signatures[[sig_name]]
  present <- intersect(genes, available_genes)

  if (length(present) >= 2) {
    tryCatch({
      seu <- AddModuleScore(seu, features = list(present), name = sig_name,
                            seed = 42, nbin = 24, ctrl = 100)
      # Rename the column (AddModuleScore adds "1" suffix)
      colnames(seu@meta.data)[ncol(seu@meta.data)] <- sig_name
      scored_sigs <- c(scored_sigs, sig_name)
      cat(sprintf("  Scored: %s (%d genes)\n", sig_name, length(present)))
    }, error = function(e) {
      cat(sprintf("  ERROR scoring %s: %s\n", sig_name, e$message))
    })
  } else {
    cat(sprintf("  SKIPPED: %s (only %d genes)\n", sig_name, length(present)))
  }
}

meta <- seu@meta.data

# Get signature columns that were scored
sig_cols <- intersect(names(signatures), colnames(meta))
cat(sprintf("\nSignatures scored: %d\n", length(sig_cols)))

# Create output directory
dir.create(OUT_TABLES_SNRNASEQ, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIGURES_SNRNASEQ, recursive = TRUE, showWarnings = FALSE)

# Summary statistics by condition
cat("\nCalculating summary statistics by condition...\n")
summary_list <- list()
for (sig in sig_cols) {
  cond_summary <- meta %>%
    group_by(Condition) %>%
    summarize(
      n_cells = n(),
      mean_score = mean(.data[[sig]], na.rm = TRUE),
      sd_score = sd(.data[[sig]], na.rm = TRUE),
      median_score = median(.data[[sig]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(signature = sig)
  summary_list[[sig]] <- cond_summary
}

summary_df <- bind_rows(summary_list)
print(summary_df)

write.csv(summary_df, file.path(OUT_TABLES_SNRNASEQ, "signature_scores_summary.csv"), row.names = FALSE)

# Statistical comparisons (Wilcoxon)
cat("\nStatistical comparisons (Implant vs Control):\n")
stats_list <- list()
for (sig in sig_cols) {
  ctrl_vals <- meta[[sig]][meta$Condition == "Control"]
  impl_vals <- meta[[sig]][meta$Condition == "Implant"]
  stab_vals <- meta[[sig]][meta$Condition == "Stab"]

  if (length(ctrl_vals) > 10 && length(impl_vals) > 10) {
    test_ic <- wilcox.test(impl_vals, ctrl_vals)
    test_sc <- wilcox.test(stab_vals, ctrl_vals)

    stats_list[[sig]] <- data.frame(
      signature = sig,
      impl_vs_ctrl_pval = test_ic$p.value,
      stab_vs_ctrl_pval = test_sc$p.value,
      impl_mean = mean(impl_vals, na.rm = TRUE),
      ctrl_mean = mean(ctrl_vals, na.rm = TRUE),
      stab_mean = mean(stab_vals, na.rm = TRUE),
      impl_vs_ctrl_fc = mean(impl_vals, na.rm = TRUE) / mean(ctrl_vals, na.rm = TRUE)
    )

    if (test_ic$p.value < 0.05) {
      cat(sprintf("  %s: p = %.2e (Implant %.3f vs Control %.3f)\n",
                  sig, test_ic$p.value, mean(impl_vals), mean(ctrl_vals)))
    }
  }
}

stats_df <- bind_rows(stats_list)
if (nrow(stats_df) > 0) {
  stats_df$impl_vs_ctrl_fdr <- p.adjust(stats_df$impl_vs_ctrl_pval, method = "BH")
  stats_df$stab_vs_ctrl_fdr <- p.adjust(stats_df$stab_vs_ctrl_pval, method = "BH")
}

write.csv(stats_df, file.path(OUT_TABLES_SNRNASEQ, "signature_statistics.csv"), row.names = FALSE)

# Dysfunction categorization based on Synaptic_Core
if ("Synaptic_Core" %in% sig_cols) {
  ctrl_mean <- mean(meta$Synaptic_Core[meta$Condition == "Control"], na.rm = TRUE)
  ctrl_sd <- sd(meta$Synaptic_Core[meta$Condition == "Control"], na.rm = TRUE)

  meta$pct_of_control <- (meta$Synaptic_Core / ctrl_mean) * 100

  meta$dysfunction_category <- case_when(
    meta$pct_of_control < 30 ~ "Catastrophic",
    meta$pct_of_control < 70 ~ "Chronic",
    meta$pct_of_control > 80 ~ "Preserved",
    TRUE ~ "Intermediate"
  )

  dys_dist <- meta %>%
    group_by(Condition, dysfunction_category) %>%
    summarize(n = n(), .groups = "drop") %>%
    group_by(Condition) %>%
    mutate(pct = 100 * n / sum(n))

  cat("\nDysfunction category distribution:\n")
  print(dys_dist)
  write.csv(dys_dist, file.path(OUT_TABLES_SNRNASEQ, "dysfunction_category_distribution.csv"), row.names = FALSE)
}

# Figures
# 1. Signature heatmap
if (length(sig_cols) > 0) {
  sig_means <- summary_df %>%
    dplyr::select(Condition, signature, mean_score) %>%
    pivot_wider(names_from = Condition, values_from = mean_score)

  sig_matrix <- as.matrix(sig_means[, -1])
  rownames(sig_matrix) <- sig_means$signature

  # Normalize to Control
  if ("Control" %in% colnames(sig_matrix)) {
    sig_norm <- sig_matrix
    for (i in 1:nrow(sig_norm)) {
      ctrl_val <- sig_matrix[i, "Control"]
      if (!is.na(ctrl_val) && abs(ctrl_val) > 0.001) {
        sig_norm[i, ] <- (sig_matrix[i, ] - ctrl_val) / abs(ctrl_val)
      }
    }

    # Heatmap as ggplot
    heatmap_df <- as.data.frame(sig_norm) %>%
      mutate(signature = rownames(sig_norm)) %>%
      pivot_longer(-signature, names_to = "Condition", values_to = "normalized_score")

    heatmap_df$Condition <- factor(heatmap_df$Condition, levels = c("Control", "Implant", "Stab"))

    p1 <- ggplot(heatmap_df, aes(x = Condition, y = signature, fill = normalized_score)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.2f", normalized_score)), size = 3) +
      scale_fill_gradient2(low = COL_DOWN, mid = "white", high = COL_UP,
                           midpoint = 0, name = "Normalized\nScore") +
      labs(title = "Neuronal Signature Scores",
           subtitle = "Normalized to Control mean",
           x = NULL, y = NULL) +
      theme_publication()

    save_figure(file.path(OUT_FIGURES_SNRNASEQ, "signature_heatmap.png"), p1, width = 7, height = 8)
  }
}

# 2. Bar plot of significant signatures
if (nrow(stats_df) > 0) {
  sig_stats <- stats_df[stats_df$impl_vs_ctrl_fdr < 0.05, ]
  if (nrow(sig_stats) > 0) {
    sig_stats$direction <- ifelse(sig_stats$impl_vs_ctrl_fc > 1, "Up", "Down")

    p2 <- ggplot(sig_stats, aes(x = reorder(signature, impl_vs_ctrl_fc), y = impl_vs_ctrl_fc - 1, fill = direction)) +
      geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.4) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      coord_flip() +
      scale_fill_manual(values = c("Up" = COL_UP, "Down" = COL_DOWN)) +
      labs(title = "Significant Signature Changes",
           subtitle = "Implant vs Control (FDR < 0.05)",
           x = NULL, y = "Fold Change - 1") +
      theme_publication() +
      theme(legend.position = "none")

    save_figure(file.path(OUT_FIGURES_SNRNASEQ, "signature_significant.png"), p2, width = 7, height = 6)
  }
}

# 3. Violin plots by condition
if (length(sig_cols) >= 3) {
  plot_sigs <- sig_cols[1:min(6, length(sig_cols))]

  violin_df <- meta %>%
    dplyr::select(Condition, all_of(plot_sigs)) %>%
    pivot_longer(-Condition, names_to = "signature", values_to = "score")

  violin_df$Condition <- factor(violin_df$Condition, levels = c("Control", "Implant", "Stab"))

  p3 <- ggplot(violin_df, aes(x = Condition, y = score, fill = Condition)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    facet_wrap(~signature, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("Control" = COL_CONTROL, "Implant" = COL_IMPLANT, "Stab" = COL_STAB)) +
    labs(title = "Signature Score Distributions by Condition",
         x = NULL, y = "Module Score") +
    theme_publication() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

  save_figure(file.path(OUT_FIGURES_SNRNASEQ, "signature_violins.png"), p3, width = 10, height = 8)
}

cat(sprintf("\nSaved signature results to: %s/\n", OUT_TABLES_SNRNASEQ))
cat(sprintf("Saved signature figures to: %s/\n", OUT_FIGURES_SNRNASEQ))
