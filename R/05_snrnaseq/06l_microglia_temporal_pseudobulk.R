# 06l_microglia_temporal_pseudobulk.R
# Temporal dynamics of SPP1+ vs Complement+ microglial populations
# Uses PSEUDOBULK (sample-level) statistics throughout
# Key question: Do SPP1+ peak early and decline (matching Jariwala rod-like)?

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/config.R")

cat("=== Microglial Population Temporal Dynamics (Pseudobulk) ===\n\n")

dir.create(file.path(OUT_TABLES_SNRNASEQ, "microglia_temporal"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_FIGURES_SNRNASEQ, "microglia_temporal"), recursive = TRUE, showWarnings = FALSE)

# Load snRNA-seq data from original source (full gene set)
SNRNASEQ_SOURCE <- SNRNASEQ_PATH  # defined in config.R
cat("Loading snRNA-seq data from source...\n")
seu <- readRDS(SNRNASEQ_SOURCE)
DefaultAssay(seu) <- "RNA"

cat(sprintf("Total cells: %d\n", ncol(seu)))
cat(sprintf("Total genes: %d\n", nrow(seu)))

# Check metadata
cat("\nMetadata columns:\n")
print(colnames(seu@meta.data))

# Identify condition and timepoint columns
meta <- seu@meta.data
cat("\nCondition distribution:\n")
print(table(meta$Condition))

# Check for Duration/Timepoint field
if ("Duration" %in% colnames(meta)) {
  cat("\nDuration distribution:\n")
  print(table(meta$Duration))
  timepoint_col <- "Duration"
} else if ("Timepoint" %in% colnames(meta)) {
  cat("\nTimepoint distribution:\n")
  print(table(meta$Timepoint))
  timepoint_col <- "Timepoint"
} else {
  cat("\nNo Duration/Timepoint column found. Checking for time info in sample names...\n")
  print(table(meta$orig.ident))
  timepoint_col <- NULL
}

# Subset to microglia
celltype_col <- if ("celltype_l2" %in% colnames(meta)) "celltype_l2" else "celltype_l1"
cat(sprintf("\nUsing celltype column: %s\n", celltype_col))
cat("Cell types present:\n")
print(table(meta[[celltype_col]]))

microglia_types <- unique(meta[[celltype_col]])[grepl("Microglia|microglia|MG", unique(meta[[celltype_col]]), ignore.case = TRUE)]
cat(sprintf("\nMicroglia types: %s\n", paste(microglia_types, collapse = ", ")))

seu_mg <- subset(seu, cells = colnames(seu)[meta[[celltype_col]] %in% microglia_types])
cat(sprintf("Microglia subset: %d cells\n", ncol(seu_mg)))

# ============================================================
# Define SPP1+ and Complement+ populations
# ============================================================

# Get gene expression
available_genes <- rownames(seu_mg)

# SPP1 module genes
spp1_genes <- c("Spp1", "Gpnmb", "Lgals3", "Fabp5", "Igf1", "Cd63", "Lpl")
spp1_present <- intersect(tools::toTitleCase(tolower(spp1_genes)), available_genes)
cat(sprintf("\nSPP1 module genes present: %s\n", paste(spp1_present, collapse = ", ")))

# Complement genes
comp_genes <- c("C1qa", "C1qb", "C1qc", "C3")
comp_present <- intersect(tools::toTitleCase(tolower(comp_genes)), available_genes)
cat(sprintf("Complement genes present: %s\n", paste(comp_present, collapse = ", ")))

# Score modules
if (length(spp1_present) >= 3) {
  seu_mg <- AddModuleScore(seu_mg, features = list(spp1_present), name = "SPP1_Module", seed = 42)
  colnames(seu_mg@meta.data)[ncol(seu_mg@meta.data)] <- "SPP1_Module"
}

if (length(comp_present) >= 3) {
  seu_mg <- AddModuleScore(seu_mg, features = list(comp_present), name = "Complement", seed = 42)
  colnames(seu_mg@meta.data)[ncol(seu_mg@meta.data)] <- "Complement"
}

# Get SPP1 expression directly
spp1_expr <- GetAssayData(seu_mg, layer = "data")["Spp1", ]
seu_mg$spp1_pos <- spp1_expr > 0

# Define Complement+ (top 25%)
comp_q75 <- quantile(seu_mg$Complement, 0.75)
seu_mg$comp_high <- seu_mg$Complement >= comp_q75

# Create population labels
seu_mg$population <- case_when(
  seu_mg$spp1_pos & seu_mg$comp_high ~ "SPP1+/Comp+",
  seu_mg$spp1_pos & !seu_mg$comp_high ~ "SPP1+/Comp-",
  !seu_mg$spp1_pos & seu_mg$comp_high ~ "SPP1-/Comp+",
  TRUE ~ "SPP1-/Comp-"
)

cat("\nOverall population distribution:\n")
print(table(seu_mg$population))

# ============================================================
# PSEUDOBULK: Aggregate by sample
# ============================================================

mg_meta <- seu_mg@meta.data

# Create sample-level summary
sample_summary <- mg_meta %>%
  group_by(orig.ident, Condition) %>%
  summarize(
    n_microglia = n(),
    n_spp1_pos = sum(spp1_pos),
    n_comp_high = sum(comp_high),
    n_spp1_only = sum(population == "SPP1+/Comp-"),
    n_comp_only = sum(population == "SPP1-/Comp+"),
    n_both = sum(population == "SPP1+/Comp+"),
    n_neither = sum(population == "SPP1-/Comp-"),
    mean_spp1_module = mean(SPP1_Module, na.rm = TRUE),
    mean_complement = mean(Complement, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_spp1_pos = 100 * n_spp1_pos / n_microglia,
    pct_comp_high = 100 * n_comp_high / n_microglia,
    pct_spp1_only = 100 * n_spp1_only / n_microglia,
    pct_comp_only = 100 * n_comp_only / n_microglia
  )

# Add timepoint if available
if (!is.null(timepoint_col)) {
  timepoint_map <- mg_meta %>%
    select(orig.ident, all_of(timepoint_col)) %>%
    distinct()
  sample_summary <- sample_summary %>%
    left_join(timepoint_map, by = "orig.ident")
  names(sample_summary)[names(sample_summary) == timepoint_col] <- "Timepoint"
}

cat("\n=== PSEUDOBULK SAMPLE SUMMARY ===\n")
print(sample_summary)

write.csv(sample_summary,
          file.path(OUT_TABLES_SNRNASEQ, "microglia_temporal", "pseudobulk_sample_summary.csv"),
          row.names = FALSE)

# ============================================================
# Statistical comparisons (sample-level)
# ============================================================

cat("\n=== SAMPLE-LEVEL STATISTICS ===\n")

# SPP1+ percentage by condition
cat("\n[1] SPP1+ percentage by Condition (sample-level):\n")
spp1_by_cond <- sample_summary %>%
  group_by(Condition) %>%
  summarize(
    n_samples = n(),
    mean_pct_spp1 = mean(pct_spp1_pos),
    sd_pct_spp1 = sd(pct_spp1_pos),
    se_pct_spp1 = sd(pct_spp1_pos) / sqrt(n()),
    .groups = "drop"
  )
print(spp1_by_cond)

# Statistical test: Implant vs Control
ctrl_spp1 <- sample_summary$pct_spp1_pos[sample_summary$Condition == "Control"]
impl_spp1 <- sample_summary$pct_spp1_pos[sample_summary$Condition == "Implant"]
stab_spp1 <- sample_summary$pct_spp1_pos[sample_summary$Condition == "Stab"]

if (length(ctrl_spp1) >= 2 && length(impl_spp1) >= 2) {
  test_impl <- wilcox.test(impl_spp1, ctrl_spp1)
  cat(sprintf("\nImplant vs Control (SPP1+): p = %.4f (Wilcoxon)\n", test_impl$p.value))
  cat(sprintf("  Control: %.2f%% ± %.2f%% (n=%d samples)\n",
              mean(ctrl_spp1), sd(ctrl_spp1), length(ctrl_spp1)))
  cat(sprintf("  Implant: %.2f%% ± %.2f%% (n=%d samples)\n",
              mean(impl_spp1), sd(impl_spp1), length(impl_spp1)))
  cat(sprintf("  Fold enrichment: %.1fx\n", mean(impl_spp1) / mean(ctrl_spp1)))
}

if (length(ctrl_spp1) >= 2 && length(stab_spp1) >= 2) {
  test_stab <- wilcox.test(stab_spp1, ctrl_spp1)
  cat(sprintf("\nStab vs Control (SPP1+): p = %.4f (Wilcoxon)\n", test_stab$p.value))
  cat(sprintf("  Stab: %.2f%% ± %.2f%% (n=%d samples)\n",
              mean(stab_spp1), sd(stab_spp1), length(stab_spp1)))
}

# Complement+ percentage by condition
cat("\n[2] Complement+ percentage by Condition (sample-level):\n")
comp_by_cond <- sample_summary %>%
  group_by(Condition) %>%
  summarize(
    n_samples = n(),
    mean_pct_comp = mean(pct_comp_high),
    sd_pct_comp = sd(pct_comp_high),
    se_pct_comp = sd(pct_comp_high) / sqrt(n()),
    .groups = "drop"
  )
print(comp_by_cond)

ctrl_comp <- sample_summary$pct_comp_high[sample_summary$Condition == "Control"]
impl_comp <- sample_summary$pct_comp_high[sample_summary$Condition == "Implant"]

if (length(ctrl_comp) >= 2 && length(impl_comp) >= 2) {
  test_comp <- wilcox.test(impl_comp, ctrl_comp)
  cat(sprintf("\nImplant vs Control (Complement+): p = %.4f (Wilcoxon)\n", test_comp$p.value))
}

# ============================================================
# Temporal analysis (if timepoint available)
# ============================================================

if ("Timepoint" %in% colnames(sample_summary)) {
  cat("\n=== TEMPORAL DYNAMICS (PSEUDOBULK) ===\n")

  # Order timepoints
  timepoint_order <- c("0wk", "1wk", "4wk", "8wk", "12wk")
  available_timepoints <- intersect(timepoint_order, unique(sample_summary$Timepoint))
  sample_summary$Timepoint <- factor(sample_summary$Timepoint, levels = available_timepoints)

  # Summary by condition and timepoint
  temporal_summary <- sample_summary %>%
    group_by(Condition, Timepoint) %>%
    summarize(
      n_samples = n(),
      mean_pct_spp1 = mean(pct_spp1_pos),
      se_pct_spp1 = sd(pct_spp1_pos) / sqrt(n()),
      mean_pct_comp = mean(pct_comp_high),
      se_pct_comp = sd(pct_comp_high) / sqrt(n()),
      mean_spp1_module = mean(mean_spp1_module),
      mean_complement_module = mean(mean_complement),
      .groups = "drop"
    )

  cat("\nTemporal summary:\n")
  print(temporal_summary)

  write.csv(temporal_summary,
            file.path(OUT_TABLES_SNRNASEQ, "microglia_temporal", "temporal_summary_pseudobulk.csv"),
            row.names = FALSE)

  # Test for temporal trend in Implant condition
  impl_temporal <- sample_summary %>% filter(Condition == "Implant")
  if (nrow(impl_temporal) >= 4) {
    # Convert timepoint to numeric for correlation
    impl_temporal$time_numeric <- as.numeric(gsub("wk", "", impl_temporal$Timepoint))

    cor_spp1 <- cor.test(impl_temporal$time_numeric, impl_temporal$pct_spp1_pos, method = "spearman")
    cor_comp <- cor.test(impl_temporal$time_numeric, impl_temporal$pct_comp_high, method = "spearman")

    cat(sprintf("\nTemporal correlation (Implant, Spearman):\n"))
    cat(sprintf("  SPP1+ vs time: rho = %.3f, p = %.4f\n", cor_spp1$estimate, cor_spp1$p.value))
    cat(sprintf("  Complement+ vs time: rho = %.3f, p = %.4f\n", cor_comp$estimate, cor_comp$p.value))
  }

  # Visualizations
  # SPP1+ over time by condition
  p_spp1_temporal <- ggplot(temporal_summary,
                            aes(x = Timepoint, y = mean_pct_spp1,
                                color = Condition, group = Condition)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_pct_spp1 - se_pct_spp1,
                      ymax = mean_pct_spp1 + se_pct_spp1), width = 0.2) +
    scale_color_manual(values = c("Control" = "#4DAF4A", "Implant" = "#E41A1C", "Stab" = "#377EB8")) +
    labs(title = "SPP1+ Microglia Over Time (Pseudobulk)",
         subtitle = "Sample-level means ± SE",
         x = "Timepoint", y = "% SPP1+ Microglia") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(OUT_FIGURES_SNRNASEQ, "microglia_temporal", "spp1_temporal_pseudobulk.png"),
         p_spp1_temporal, width = 8, height = 5)

  # Complement+ over time
  p_comp_temporal <- ggplot(temporal_summary,
                            aes(x = Timepoint, y = mean_pct_comp,
                                color = Condition, group = Condition)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_pct_comp - se_pct_comp,
                      ymax = mean_pct_comp + se_pct_comp), width = 0.2) +
    scale_color_manual(values = c("Control" = "#4DAF4A", "Implant" = "#E41A1C", "Stab" = "#377EB8")) +
    labs(title = "Complement+ Microglia Over Time (Pseudobulk)",
         subtitle = "Sample-level means ± SE",
         x = "Timepoint", y = "% Complement-High Microglia") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(OUT_FIGURES_SNRNASEQ, "microglia_temporal", "complement_temporal_pseudobulk.png"),
         p_comp_temporal, width = 8, height = 5)

  # Combined panel
  temporal_long <- temporal_summary %>%
    select(Condition, Timepoint, mean_pct_spp1, mean_pct_comp) %>%
    pivot_longer(cols = c(mean_pct_spp1, mean_pct_comp),
                 names_to = "Population", values_to = "Percentage") %>%
    mutate(Population = recode(Population,
                               "mean_pct_spp1" = "SPP1+ (rod-like)",
                               "mean_pct_comp" = "Complement+"))

  p_combined <- ggplot(temporal_long,
                       aes(x = Timepoint, y = Percentage,
                           color = Population, linetype = Condition, group = interaction(Population, Condition))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("SPP1+ (rod-like)" = "#E41A1C", "Complement+" = "#984EA3")) +
    facet_wrap(~Condition, ncol = 3) +
    labs(title = "Microglial Population Dynamics",
         subtitle = "SPP1+/rod-like vs Complement+ over time (pseudobulk)",
         x = "Timepoint", y = "% of Microglia") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(OUT_FIGURES_SNRNASEQ, "microglia_temporal", "combined_temporal_pseudobulk.png"),
         p_combined, width = 12, height = 5)
}

# ============================================================
# SUMMARY
# ============================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("MICROGLIAL TEMPORAL DYNAMICS SUMMARY (PSEUDOBULK)\n")
cat(strrep("=", 70), "\n")

cat("\n[1] SPP1+ MICROGLIA (rod-like correlate):\n")
cat(sprintf("    Control: %.2f%% ± %.2f%%\n", mean(ctrl_spp1), sd(ctrl_spp1)))
cat(sprintf("    Implant: %.2f%% ± %.2f%% (%.1fx enrichment)\n",
            mean(impl_spp1), sd(impl_spp1), mean(impl_spp1)/mean(ctrl_spp1)))
if (exists("test_impl")) {
  cat(sprintf("    Implant vs Control: p = %.4f\n", test_impl$p.value))
}

cat("\n[2] COMPLEMENT+ MICROGLIA:\n")
cat(sprintf("    Control: %.2f%% ± %.2f%%\n", mean(ctrl_comp), sd(ctrl_comp)))
cat(sprintf("    Implant: %.2f%% ± %.2f%%\n", mean(impl_comp), sd(impl_comp)))
if (exists("test_comp")) {
  cat(sprintf("    Implant vs Control: p = %.4f\n", test_comp$p.value))
}

cat("\n[3] KEY QUESTION - Does SPP1+ peak early and decline?\n")
cat("    (This would match Jariwala rod-like morphology finding)\n")
if (exists("cor_spp1")) {
  if (cor_spp1$estimate < -0.3 && cor_spp1$p.value < 0.1) {
    cat("    RESULT: SPP1+ shows DECLINING trend over time (supports Jariwala)\n")
  } else if (cor_spp1$estimate > 0.3 && cor_spp1$p.value < 0.1) {
    cat("    RESULT: SPP1+ shows INCREASING trend over time (contradicts Jariwala)\n")
  } else {
    cat("    RESULT: No significant temporal trend detected\n")
  }
  cat(sprintf("    Spearman rho = %.3f, p = %.4f\n", cor_spp1$estimate, cor_spp1$p.value))
}

cat("\n[4] COMPLEMENT+ TEMPORAL TREND:\n")
if (exists("cor_comp")) {
  if (cor_comp$estimate > 0.3 && cor_comp$p.value < 0.1) {
    cat("    RESULT: Complement+ INCREASES over time (maladaptive persistence)\n")
  } else if (cor_comp$estimate < -0.3 && cor_comp$p.value < 0.1) {
    cat("    RESULT: Complement+ DECREASES over time\n")
  } else {
    cat("    RESULT: No significant temporal trend detected\n")
  }
  cat(sprintf("    Spearman rho = %.3f, p = %.4f\n", cor_comp$estimate, cor_comp$p.value))
}

cat("\nOutputs saved to:\n")
cat(sprintf("  %s/microglia_temporal/\n", OUT_TABLES_SNRNASEQ))
cat(sprintf("  %s/microglia_temporal/\n", OUT_FIGURES_SNRNASEQ))
