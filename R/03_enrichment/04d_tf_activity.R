# 04d_tf_activity.R
# Transcription Factor Activity Inference using dorothea + viper
# Input: DEG results
# Output: tables/enrichment/tf_*.csv, figures/enrichment/tf_*.png

library(dorothea)
library(viper)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/config.R")

# Create output directories
dir.create(file.path(OUT_TABLES, "enrichment"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_FIGURES, "enrichment"), recursive = TRUE, showWarnings = FALSE)

# Load DEG results
cat("Loading DEG results...\n")
silicon <- read.csv(file.path(OUT_TABLES_DEG, "silicon_deseq2_results.csv"))
polyimide <- read.csv(file.path(OUT_TABLES_DEG, "polyimide_limma_results.csv"))

# Get DoRothEA regulons (mouse, convert to rat)
cat("Loading DoRothEA regulons...\n")
regulons_raw <- dorothea_mm  # Mouse regulons from dorothea package
regulons_filtered <- regulons_raw %>%
  filter(confidence %in% c("A", "B", "C"))

# Convert to title case for rat (mouse genes are Title case, same as rat)
cat(sprintf("Loaded regulons: %d TF-target pairs\n", nrow(regulons_filtered)))
cat(sprintf("TFs: %d, Targets: %d\n",
            length(unique(regulons_filtered$tf)),
            length(unique(regulons_filtered$target))))

# Convert to viper regulon format
regulon_list <- split(regulons_filtered, regulons_filtered$tf)
viper_regulon <- lapply(regulon_list, function(df) {
  list(
    tfmode = setNames(df$mor, df$target),
    likelihood = rep(1, nrow(df))
  )
})
class(viper_regulon) <- "regulon"

# Prepare expression matrices
cat("\nPreparing expression matrices...\n")

# Silicon
silicon_mat <- silicon %>%
  filter(!is.na(log2FoldChange) & !is.na(gene) & gene != "") %>%
  group_by(gene) %>%
  slice_max(abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup()

silicon_expr <- setNames(silicon_mat$log2FoldChange, silicon_mat$gene)
silicon_expr_mat <- matrix(silicon_expr, ncol = 1, dimnames = list(names(silicon_expr), "Silicon"))

# Polyimide
polyimide_mat <- polyimide %>%
  filter(!is.na(logFC) & !is.na(gene) & gene != "") %>%
  group_by(gene) %>%
  slice_max(abs(logFC), n = 1, with_ties = FALSE) %>%
  ungroup()

polyimide_expr <- setNames(polyimide_mat$logFC, polyimide_mat$gene)
polyimide_expr_mat <- matrix(polyimide_expr, ncol = 1, dimnames = list(names(polyimide_expr), "Polyimide"))

# Run viper for TF activity
cat("\nRunning TF activity inference for Silicon...\n")
tf_silicon <- viper(silicon_expr_mat, viper_regulon, minsize = 5, eset.filter = FALSE, verbose = FALSE)
tf_silicon_df <- data.frame(
  source = rownames(tf_silicon),
  score = as.numeric(tf_silicon[, 1]),
  platform = "Silicon"
) %>%
  arrange(desc(abs(score)))

cat(sprintf("Silicon: %d TFs with activity scores\n", nrow(tf_silicon_df)))

cat("Running TF activity inference for Polyimide...\n")
tf_polyimide <- viper(polyimide_expr_mat, viper_regulon, minsize = 5, eset.filter = FALSE, verbose = FALSE)
tf_polyimide_df <- data.frame(
  source = rownames(tf_polyimide),
  score = as.numeric(tf_polyimide[, 1]),
  platform = "Polyimide"
) %>%
  arrange(desc(abs(score)))

cat(sprintf("Polyimide: %d TFs with activity scores\n", nrow(tf_polyimide_df)))

# Save results
write.csv(tf_silicon_df, file.path(OUT_TABLES, "enrichment", "tf_activity_silicon.csv"), row.names = FALSE)
write.csv(tf_polyimide_df, file.path(OUT_TABLES, "enrichment", "tf_activity_polyimide.csv"), row.names = FALSE)

# Cross-platform TF concordance
common_tfs <- intersect(tf_silicon_df$source, tf_polyimide_df$source)
cat(sprintf("\nCommon TFs: %d\n", length(common_tfs)))

if (length(common_tfs) >= 5) {
  tf_merged <- merge(
    tf_silicon_df[, c("source", "score")],
    tf_polyimide_df[, c("source", "score")],
    by = "source", suffixes = c("_silicon", "_polyimide")
  )

  cor_test <- cor.test(tf_merged$score_silicon, tf_merged$score_polyimide, method = "spearman")
  cat(sprintf("TF activity correlation: rho = %.3f (p = %.2e)\n", cor_test$estimate, cor_test$p.value))

  # TF concordance plot
  p1 <- ggplot(tf_merged, aes(x = score_silicon, y = score_polyimide)) +
    geom_point(aes(color = score_silicon * score_polyimide > 0), alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
    geom_hline(yintercept = 0, color = COL_REF, linewidth = 0.5) +
    geom_vline(xintercept = 0, color = COL_REF, linewidth = 0.5) +
    ggrepel::geom_text_repel(
      data = tf_merged %>% filter(abs(score_silicon) > 2 | abs(score_polyimide) > 2),
      aes(label = source), size = 3, max.overlaps = 15
    ) +
    scale_color_manual(values = c("FALSE" = COL_NS, "TRUE" = COL_SIG),
                       labels = c("Discordant", "Concordant")) +
    labs(title = "Transcription Factor Activity Concordance",
         subtitle = sprintf("Spearman rho = %.2f, n = %d TFs", cor_test$estimate, nrow(tf_merged)),
         x = "TF Activity (Silicon)", y = "TF Activity (Polyimide)", color = "Direction") +
    theme_publication()

  save_figure(file.path(OUT_FIGURES, "enrichment", "tf_activity_concordance.png"), p1, width = 8, height = 7)

  # Save concordance stats
  # NOTE: TF activity score comparison is valid across platforms (normalized scores)
  tf_concordance <- data.frame(
    metric = "TF_activity_correlation",
    rank_correlation = cor_test$estimate,
    p_value = cor_test$p.value,
    n_tfs = length(common_tfs)
  )
  write.csv(tf_concordance, file.path(OUT_TABLES, "enrichment", "tf_concordance.csv"), row.names = FALSE)
}

# Top TFs bar plot
top_silicon <- tf_silicon_df %>%
  arrange(desc(abs(score))) %>%
  head(20)

top_polyimide <- tf_polyimide_df %>%
  arrange(desc(abs(score))) %>%
  head(20)

top_tfs <- bind_rows(top_silicon, top_polyimide) %>%
  mutate(direction = ifelse(score > 0, "Activated", "Repressed"))

p2 <- ggplot(top_tfs, aes(x = reorder(source, score), y = score, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  facet_wrap(~platform, scales = "free_y") +
  coord_flip() +
  scale_fill_manual(values = c("Activated" = COL_UP, "Repressed" = COL_DOWN)) +
  labs(title = "Top Transcription Factor Activities",
       subtitle = "DoRothEA regulons (confidence A-C)",
       x = NULL, y = "Normalized Activity Score") +
  theme_publication() +
  theme(legend.position = "bottom")

save_figure(file.path(OUT_FIGURES, "enrichment", "tf_top_activities.png"), p2, width = 12, height = 10)

# Key inflammatory and neuronal TFs
inflammatory_tfs <- c("Nfkb1", "Rela", "Stat1", "Stat3", "Irf1", "Irf7", "Jun", "Fos")
neuronal_tfs <- c("Creb1", "Mef2a", "Mef2c", "Neurod1", "Rest", "Nfia")

key_tfs <- tf_merged %>%
  filter(source %in% c(inflammatory_tfs, neuronal_tfs)) %>%
  mutate(category = ifelse(source %in% inflammatory_tfs, "Inflammatory", "Neuronal"))

if (nrow(key_tfs) > 0) {
  key_tfs_long <- key_tfs %>%
    pivot_longer(cols = c(score_silicon, score_polyimide),
                 names_to = "platform", values_to = "activity") %>%
    mutate(platform = gsub("score_", "", platform))

  p3 <- ggplot(key_tfs_long, aes(x = source, y = activity, fill = platform)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    facet_wrap(~category, scales = "free_x") +
    scale_fill_manual(values = c("silicon" = COL_IMPLANT, "polyimide" = COL_STAB)) +
    labs(title = "Key Transcription Factor Activities",
         x = NULL, y = "Activity Score", fill = "Platform") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_figure(file.path(OUT_FIGURES, "enrichment", "tf_key_activities.png"), p3, width = 10, height = 6)
}

cat(sprintf("\nSaved TF activity results to: %s/enrichment/\n", OUT_TABLES))
