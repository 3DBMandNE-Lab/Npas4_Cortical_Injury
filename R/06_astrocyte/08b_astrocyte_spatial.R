# 08b_astrocyte_spatial.R
# Spatial analysis of astrocyte signatures and C3 expression
# Tests the two-cell complement model: Microglia (C1q) vs Astrocytes (C3)

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
source("R/config.R")

cat("=== Astrocyte Spatial Analysis: Two-Cell Complement Model ===\n\n")

dir.create(file.path(OUT_FIGURES_MANUSCRIPT, "astrocyte_spatial"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_TABLES_SPATIAL, "astrocyte"), recursive = TRUE, showWarnings = FALSE)

# Load spatial data
SPATIAL_DIR <- "data/external/stRNAseq"
sample_names <- c("Visium_4A", "Visium_4B", "Visium_5B", "Visium_7C", "Visium_8A", "Visium_8C")

cat("Loading Visium data from source files...\n")
spatial_list <- list()

if (TRUE) {  # Always load from source

  for (sname in sample_names) {
    sample_dir <- file.path(SPATIAL_DIR, sname)
    matrix_dir <- file.path(sample_dir, "filtered_feature_bc_matrix")
    spatial_dir <- file.path(sample_dir, "spatial")

    if (dir.exists(matrix_dir)) {
      cat(sprintf("  Loading %s...\n", sname))
      counts <- Read10X(matrix_dir)
      seu <- CreateSeuratObject(counts = counts, project = sname, assay = "Spatial")

      positions_file <- file.path(spatial_dir, "tissue_positions.csv")
      if (file.exists(positions_file)) {
        positions <- read.csv(positions_file, header = TRUE)
        if ("barcode" %in% colnames(positions)) {
          rownames(positions) <- positions$barcode
        }
        common_cells <- intersect(colnames(seu), rownames(positions))
        positions <- positions[common_cells, ]
        seu$array_row <- positions[colnames(seu), "array_row"]
        seu$array_col <- positions[colnames(seu), "array_col"]
      }

      seu <- NormalizeData(seu, verbose = FALSE)
      seu$sample <- sname
      spatial_list[[sname]] <- seu
      cat(sprintf("    Loaded %d spots\n", ncol(seu)))
    }
  }
}

# ============================================================
# Define gene signatures for two-cell model
# ============================================================

# Astrocyte markers (general)
astrocyte_genes <- c("Gfap", "Aqp4", "S100b", "Aldh1l1", "Slc1a2", "Slc1a3")

# Astrocyte C3 production (the key finding)
astro_c3_genes <- c("C3", "C4a", "C4b")

# Microglial C1q (initiation)
micro_c1q_genes <- c("C1qa", "C1qb", "C1qc")

# Full complement cascade
complement_full <- c("C1qa", "C1qb", "C1qc", "C3", "C4a", "C4b")

# Microglial markers
microglia_genes <- c("Cx3cr1", "P2ry12", "Tmem119", "Aif1", "Cd68")

# SPP1 module (for comparison)
spp1_genes <- c("Spp1", "Gpnmb", "Lgals3", "Fabp5", "Cd63")

results_list <- list()

for (sample_name in names(spatial_list)) {
  cat(sprintf("\n=== Processing %s ===\n", sample_name))
  seu <- spatial_list[[sample_name]]
  available_genes <- rownames(seu)

  # Score each module
  # Astrocyte general signature
  astro_present <- intersect(astrocyte_genes, available_genes)
  if (length(astro_present) >= 3) {
    seu <- AddModuleScore(seu, features = list(astro_present), name = "Astrocyte", seed = 42)
    colnames(seu@meta.data)[ncol(seu@meta.data)] <- "Astrocyte"
    cat(sprintf("  Astrocyte module: %d/%d genes\n", length(astro_present), length(astrocyte_genes)))
  }

  # Microglial signature
  micro_present <- intersect(microglia_genes, available_genes)
  if (length(micro_present) >= 3) {
    seu <- AddModuleScore(seu, features = list(micro_present), name = "Microglia", seed = 42)
    colnames(seu@meta.data)[ncol(seu@meta.data)] <- "Microglia"
    cat(sprintf("  Microglia module: %d/%d genes\n", length(micro_present), length(microglia_genes)))
  }

  # C1q module (microglial initiation)
  c1q_present <- intersect(micro_c1q_genes, available_genes)
  if (length(c1q_present) >= 2) {
    seu <- AddModuleScore(seu, features = list(c1q_present), name = "C1q_Module", seed = 42)
    colnames(seu@meta.data)[ncol(seu@meta.data)] <- "C1q_Module"
    cat(sprintf("  C1q module: %d/%d genes\n", length(c1q_present), length(micro_c1q_genes)))
  }

  # C3 expression (key for astrocyte contribution)
  if ("C3" %in% available_genes) {
    seu$C3_expr <- GetAssayData(seu, layer = "data")["C3", ]
    cat(sprintf("  C3 expression: mean = %.3f\n", mean(seu$C3_expr)))
  }

  # Gfap expression (astrocyte marker)
  if ("Gfap" %in% available_genes) {
    seu$Gfap_expr <- GetAssayData(seu, layer = "data")["Gfap", ]
    cat(sprintf("  Gfap expression: mean = %.3f\n", mean(seu$Gfap_expr)))
  }

  # Individual C1q genes
  for (gene in c("C1qa", "C1qb", "C1qc")) {
    if (gene %in% available_genes) {
      seu@meta.data[[paste0(gene, "_expr")]] <- GetAssayData(seu, layer = "data")[gene, ]
    }
  }

  # SPP1 module
  spp1_present <- intersect(spp1_genes, available_genes)
  if (length(spp1_present) >= 3) {
    seu <- AddModuleScore(seu, features = list(spp1_present), name = "SPP1_Module", seed = 42)
    colnames(seu@meta.data)[ncol(seu@meta.data)] <- "SPP1_Module"
  }

  # ============================================================
  # Key correlations for two-cell model
  # ============================================================

  correlations <- list()

  # Astrocyte vs Microglia (should be distinct)
  if ("Astrocyte" %in% colnames(seu@meta.data) && "Microglia" %in% colnames(seu@meta.data)) {
    cor_astro_micro <- cor.test(seu$Astrocyte, seu$Microglia, method = "pearson")
    correlations$astro_micro <- cor_astro_micro$estimate
    cat(sprintf("  Astrocyte vs Microglia: r = %.3f\n", cor_astro_micro$estimate))
  }

  # C3 vs GFAP (astrocyte C3 production)
  if ("C3_expr" %in% colnames(seu@meta.data) && "Gfap_expr" %in% colnames(seu@meta.data)) {
    cor_c3_gfap <- cor.test(seu$C3_expr, seu$Gfap_expr, method = "pearson")
    correlations$c3_gfap <- cor_c3_gfap$estimate
    cat(sprintf("  C3 vs GFAP: r = %.3f (astrocyte C3 production)\n", cor_c3_gfap$estimate))
  }

  # C3 vs Astrocyte module
  if ("C3_expr" %in% colnames(seu@meta.data) && "Astrocyte" %in% colnames(seu@meta.data)) {
    cor_c3_astro <- cor.test(seu$C3_expr, seu$Astrocyte, method = "pearson")
    correlations$c3_astro <- cor_c3_astro$estimate
    cat(sprintf("  C3 vs Astrocyte module: r = %.3f\n", cor_c3_astro$estimate))
  }

  # C3 vs C1q (complement cascade coherence)
  if ("C3_expr" %in% colnames(seu@meta.data) && "C1q_Module" %in% colnames(seu@meta.data)) {
    cor_c3_c1q <- cor.test(seu$C3_expr, seu$C1q_Module, method = "pearson")
    correlations$c3_c1q <- cor_c3_c1q$estimate
    cat(sprintf("  C3 vs C1q: r = %.3f (complement cascade)\n", cor_c3_c1q$estimate))
  }

  # C1q vs Microglia (microglial C1q)
  if ("C1q_Module" %in% colnames(seu@meta.data) && "Microglia" %in% colnames(seu@meta.data)) {
    cor_c1q_micro <- cor.test(seu$C1q_Module, seu$Microglia, method = "pearson")
    correlations$c1q_micro <- cor_c1q_micro$estimate
    cat(sprintf("  C1q vs Microglia: r = %.3f (microglial C1q)\n", cor_c1q_micro$estimate))
  }

  # C3 vs Microglia (microglial C3 - should be weaker than astrocyte)
  if ("C3_expr" %in% colnames(seu@meta.data) && "Microglia" %in% colnames(seu@meta.data)) {
    cor_c3_micro <- cor.test(seu$C3_expr, seu$Microglia, method = "pearson")
    correlations$c3_micro <- cor_c3_micro$estimate
    cat(sprintf("  C3 vs Microglia: r = %.3f (microglial C3)\n", cor_c3_micro$estimate))
  }

  # Store results
  results_list[[sample_name]] <- data.frame(
    sample = sample_name,
    n_spots = ncol(seu),
    r_astro_micro = correlations$astro_micro %||% NA,
    r_c3_gfap = correlations$c3_gfap %||% NA,
    r_c3_astro = correlations$c3_astro %||% NA,
    r_c3_c1q = correlations$c3_c1q %||% NA,
    r_c1q_micro = correlations$c1q_micro %||% NA,
    r_c3_micro = correlations$c3_micro %||% NA
  )

  # Update spatial list
  spatial_list[[sample_name]] <- seu
}

# Combine results
results_df <- bind_rows(results_list)
cat("\n\n=== SUMMARY: Two-Cell Complement Model Correlations ===\n")
print(results_df)

write.csv(results_df, file.path(OUT_TABLES_SPATIAL, "astrocyte", "spatial_correlation_summary.csv"), row.names = FALSE)

# ============================================================
# Key test: C3-GFAP vs C3-Microglia correlation comparison
# ============================================================

cat("\n\n=== TWO-CELL MODEL TEST ===\n")
mean_c3_astro <- mean(results_df$r_c3_gfap, na.rm = TRUE)
mean_c3_micro <- mean(results_df$r_c3_micro, na.rm = TRUE)
mean_c1q_micro <- mean(results_df$r_c1q_micro, na.rm = TRUE)

cat(sprintf("\nMean C3-GFAP correlation (astrocyte C3): r = %.3f\n", mean_c3_astro))
cat(sprintf("Mean C3-Microglia correlation: r = %.3f\n", mean_c3_micro))
cat(sprintf("Mean C1q-Microglia correlation: r = %.3f\n", mean_c1q_micro))

if (mean_c3_astro > mean_c3_micro) {
  cat("\n[SUPPORTED]: C3 correlates MORE with astrocytes than microglia spatially\n")
  cat("This supports the two-cell model:\n")
  cat("  - C1q initiation by microglia (r = %.3f with microglial markers)\n", mean_c1q_micro)
  cat("  - C3 amplification by astrocytes (r = %.3f with GFAP)\n", mean_c3_astro)
} else {
  cat("\n[NOT SUPPORTED]: C3 correlates similarly or more with microglia\n")
  cat("This does not support astrocyte-dominant C3 production spatially\n")
}

# ============================================================
# Visualizations for chronic sample (Visium_7C)
# ============================================================

chronic_sample <- "Visium_7C"
if (chronic_sample %in% names(spatial_list)) {
  cat(sprintf("\n\nGenerating visualizations for %s...\n", chronic_sample))
  seu <- spatial_list[[chronic_sample]]

  plot_data <- seu@meta.data
  plot_data$array_row <- as.numeric(plot_data$array_row)
  plot_data$array_col <- as.numeric(plot_data$array_col)

  # Panel A: GFAP (astrocyte marker)
  p1 <- ggplot(plot_data, aes(x = array_col, y = -array_row, color = Gfap_expr)) +
    geom_point(size = 1.5) +
    scale_color_gradient(low = "grey90", high = "#4DAF4A") +
    labs(title = "GFAP (Astrocytes)", x = NULL, y = NULL, color = "Expr") +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    coord_fixed()

  # Panel B: C1q (microglial initiation)
  p2 <- ggplot(plot_data, aes(x = array_col, y = -array_row, color = C1q_Module)) +
    geom_point(size = 1.5) +
    scale_color_gradient(low = "grey90", high = "#984EA3") +
    labs(title = "C1q Module (Microglia)", x = NULL, y = NULL, color = "Score") +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    coord_fixed()

  # Panel C: C3 (astrocyte amplification)
  p3 <- ggplot(plot_data, aes(x = array_col, y = -array_row, color = C3_expr)) +
    geom_point(size = 1.5) +
    scale_color_gradient(low = "grey90", high = "#E41A1C") +
    labs(title = "C3 Expression", x = NULL, y = NULL, color = "Expr") +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    coord_fixed()

  # Panel D: Microglia module
  p4 <- ggplot(plot_data, aes(x = array_col, y = -array_row, color = Microglia)) +
    geom_point(size = 1.5) +
    scale_color_gradient(low = "grey90", high = "#377EB8") +
    labs(title = "Microglia Module", x = NULL, y = NULL, color = "Score") +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    coord_fixed()

  # Combine spatial maps
  spatial_combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(title = sprintf("Two-Cell Complement Model - %s", chronic_sample))

  ggsave(file.path(OUT_FIGURES_MANUSCRIPT, "astrocyte_spatial", "two_cell_model_spatial.png"),
         spatial_combined, width = 12, height = 10)

  # Scatter plots for correlation
  # C3 vs GFAP (astrocyte C3)
  r_c3_gfap <- results_df$r_c3_gfap[results_df$sample == chronic_sample]
  p5 <- ggplot(plot_data, aes(x = Gfap_expr, y = C3_expr)) +
    geom_point(alpha = 0.3, size = 1, color = "#4DAF4A") +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
    labs(title = sprintf("C3 vs GFAP (r = %.2f)", r_c3_gfap),
         subtitle = "Astrocyte C3 production",
         x = "GFAP Expression", y = "C3 Expression") +
    theme_minimal()

  # C3 vs Microglia
  r_c3_micro <- results_df$r_c3_micro[results_df$sample == chronic_sample]
  p6 <- ggplot(plot_data, aes(x = Microglia, y = C3_expr)) +
    geom_point(alpha = 0.3, size = 1, color = "#377EB8") +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
    labs(title = sprintf("C3 vs Microglia (r = %.2f)", r_c3_micro),
         subtitle = "Microglial C3?",
         x = "Microglia Module Score", y = "C3 Expression") +
    theme_minimal()

  # C1q vs Microglia
  r_c1q_micro <- results_df$r_c1q_micro[results_df$sample == chronic_sample]
  p7 <- ggplot(plot_data, aes(x = Microglia, y = C1q_Module)) +
    geom_point(alpha = 0.3, size = 1, color = "#984EA3") +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
    labs(title = sprintf("C1q vs Microglia (r = %.2f)", r_c1q_micro),
         subtitle = "Microglial C1q initiation",
         x = "Microglia Module Score", y = "C1q Module Score") +
    theme_minimal()

  # C3 vs C1q (cascade)
  r_c3_c1q <- results_df$r_c3_c1q[results_df$sample == chronic_sample]
  p8 <- ggplot(plot_data, aes(x = C1q_Module, y = C3_expr)) +
    geom_point(alpha = 0.3, size = 1, color = "#FF7F00") +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
    labs(title = sprintf("C3 vs C1q (r = %.2f)", r_c3_c1q),
         subtitle = "Complement cascade coordination",
         x = "C1q Module Score", y = "C3 Expression") +
    theme_minimal()

  scatter_combined <- (p5 | p6) / (p7 | p8)
  ggsave(file.path(OUT_FIGURES_MANUSCRIPT, "astrocyte_spatial", "two_cell_model_correlations.png"),
         scatter_combined, width = 10, height = 8)

  cat("Saved visualizations\n")
}

# ============================================================
# Forest plot of correlations across samples
# ============================================================

# Reshape for forest plot
forest_data <- results_df %>%
  select(sample, r_c3_gfap, r_c3_micro, r_c1q_micro) %>%
  pivot_longer(-sample, names_to = "comparison", values_to = "correlation") %>%
  mutate(comparison = recode(comparison,
    "r_c3_gfap" = "C3 ~ GFAP (Astrocyte C3)",
    "r_c3_micro" = "C3 ~ Microglia (Microglial C3?)",
    "r_c1q_micro" = "C1q ~ Microglia (Microglial C1q)"
  ))

p_forest <- ggplot(forest_data, aes(x = correlation, y = sample, color = comparison)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_color_manual(values = c(
    "C3 ~ GFAP (Astrocyte C3)" = "#4DAF4A",
    "C3 ~ Microglia (Microglial C3?)" = "#377EB8",
    "C1q ~ Microglia (Microglial C1q)" = "#984EA3"
  )) +
  labs(title = "Two-Cell Complement Model: Spatial Correlations",
       subtitle = "C3 correlates more with astrocytes (GFAP) than microglia across samples",
       x = "Pearson Correlation (r)", y = NULL, color = "Comparison") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_FIGURES_MANUSCRIPT, "astrocyte_spatial", "two_cell_model_forest.png"),
       p_forest, width = 10, height = 6)

# ============================================================
# Summary statistics
# ============================================================

summary_stats <- data.frame(
  comparison = c("C3 ~ GFAP (Astrocyte)", "C3 ~ Microglia", "C1q ~ Microglia"),
  mean_r = c(
    mean(results_df$r_c3_gfap, na.rm = TRUE),
    mean(results_df$r_c3_micro, na.rm = TRUE),
    mean(results_df$r_c1q_micro, na.rm = TRUE)
  ),
  sd_r = c(
    sd(results_df$r_c3_gfap, na.rm = TRUE),
    sd(results_df$r_c3_micro, na.rm = TRUE),
    sd(results_df$r_c1q_micro, na.rm = TRUE)
  ),
  min_r = c(
    min(results_df$r_c3_gfap, na.rm = TRUE),
    min(results_df$r_c3_micro, na.rm = TRUE),
    min(results_df$r_c1q_micro, na.rm = TRUE)
  ),
  max_r = c(
    max(results_df$r_c3_gfap, na.rm = TRUE),
    max(results_df$r_c3_micro, na.rm = TRUE),
    max(results_df$r_c1q_micro, na.rm = TRUE)
  )
)

cat("\n\n=== Summary Statistics ===\n")
print(summary_stats)

write.csv(summary_stats, file.path(OUT_TABLES_SPATIAL, "astrocyte", "two_cell_model_summary.csv"), row.names = FALSE)

cat("\n\n=== INTERPRETATION ===\n")
cat("Two-cell complement model test:\n")
cat("1. C1q (initiation) should correlate with microglia - TESTED\n")
cat("2. C3 (amplification) should correlate MORE with astrocytes than microglia - TESTED\n")
cat("3. C3 and C1q should correlate (cascade coherence) - TESTED\n\n")

cat("If C3~GFAP > C3~Microglia, this SUPPORTS the two-cell model from snRNA-seq.\n")
cat("If C3~GFAP ≈ C3~Microglia, the spatial evidence is AMBIGUOUS.\n")
cat("If C3~GFAP < C3~Microglia, this CONTRADICTS the two-cell model.\n")

cat("\n\nOutputs saved to:\n")
cat(sprintf("  %s/astrocyte/\n", OUT_TABLES_SPATIAL))
cat(sprintf("  %s/astrocyte_spatial/\n", OUT_FIGURES_MANUSCRIPT))
