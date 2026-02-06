# 08_astrocyte_analysis.R
# Comprehensive astrocyte analysis to address asymmetric focus on microglia
#
# Key questions:
# 1. Reactive vs homeostatic astrocyte states by condition
# 2. Complement (C3) production - astrocytes as complement source
# 3. Glial scar markers (GFAP, Vim, Tnc)
# 4. Astrocyte-specific DE between conditions
# 5. Comparison to microglial findings
# 6. Spatial distribution (Visium)

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# Source config
source("R/config.R")

# Output directories
dir.create("output/tables/snrnaseq/astrocytes", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures/snrnaseq/astrocytes", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# PART 1: Load data and subset astrocytes
# =============================================================================

cat("Loading Seurat object...\n")
seu <- readRDS(SNRNASEQ_PATH)

# Subset to astrocytes
astro <- subset(seu, celltype_l1 == "Astrocyte")
cat("Total astrocytes:", ncol(astro), "\n")
cat("By condition:\n")
print(table(astro$Condition))

# Set default assay to RNA for module scoring
DefaultAssay(astro) <- "RNA"
cat("Using assay:", DefaultAssay(astro), "\n")

# Check available genes
cat("Genes in RNA assay:", nrow(astro), "\n")

# =============================================================================
# PART 2: Astrocyte state analysis
# =============================================================================

cat("\n=== Astrocyte State Analysis ===\n")

# Define reactive astrocyte markers (A1/neurotoxic)
reactive_markers <- c("Gfap", "Vim", "Lcn2", "Serpina3n", "C3", "Gbp2", "Psmb8",
                      "Srgn", "H2-D1", "H2-K1", "Serping1", "Fkbp5")

# Homeostatic astrocyte markers
homeostatic_markers <- c("Slc1a2", "Slc1a3", "Aldh1l1", "Aqp4", "Gja1", "Sox9",
                         "Ndrg2", "S100b", "Glul", "Gjb6")

# A2/neuroprotective markers (removed Cd14 - not in object)
neuroprotective_markers <- c("Emp1", "S100a10", "Tm4sf1", "Sphk1", "Cd109",
                             "Ptx3", "Clcf1")

# Glial scar markers
scar_markers <- c("Gfap", "Vim", "Tnc", "Cspg4", "Ncan", "Bcan")

# Complement-related (astrocytes as complement source)
complement_markers <- c("C3", "C4a", "C4b", "Serping1", "Cfh", "Cfi", "C1s", "C1r")

# Filter markers to those present in the object
filter_markers <- function(markers, obj) {
  present <- markers[markers %in% rownames(obj)]
  cat("  Using", length(present), "of", length(markers), "markers\n")
  return(present)
}

# Calculate module scores
cat("Calculating module scores...\n")
# Use simple mean expression as module score (more robust for small populations)
calc_module_score <- function(obj, markers, score_name) {
  present <- markers[markers %in% rownames(obj)]
  cat("  ", score_name, ": Using", length(present), "of", length(markers), "markers\n")
  if(length(present) < 2) return(rep(NA, ncol(obj)))

  # Get normalized expression and calculate mean
  expr <- GetAssayData(obj, layer = "data")[present, , drop = FALSE]
  return(colMeans(as.matrix(expr)))
}

astro$Astrocyte_reactive_score <- calc_module_score(astro, reactive_markers, "Reactive")
astro$Astrocyte_homeostatic_score <- calc_module_score(astro, homeostatic_markers, "Homeostatic")
astro$A2_neuroprotective_score <- calc_module_score(astro, neuroprotective_markers, "A2_neuroprotective")
astro$GlialScar_score <- calc_module_score(astro, scar_markers, "GlialScar")
astro$Complement_production_score <- calc_module_score(astro, complement_markers, "Complement")

# Summarize by condition
state_summary <- astro@meta.data %>%
  group_by(Condition) %>%
  summarise(
    n_cells = n(),
    reactive_mean = mean(Astrocyte_reactive_score, na.rm = TRUE),
    reactive_sd = sd(Astrocyte_reactive_score, na.rm = TRUE),
    homeostatic_mean = mean(Astrocyte_homeostatic_score, na.rm = TRUE),
    homeostatic_sd = sd(Astrocyte_homeostatic_score, na.rm = TRUE),
    A2_mean = mean(A2_neuroprotective_score, na.rm = TRUE),
    A2_sd = sd(A2_neuroprotective_score, na.rm = TRUE),
    scar_mean = mean(GlialScar_score, na.rm = TRUE),
    scar_sd = sd(GlialScar_score, na.rm = TRUE),
    complement_mean = mean(Complement_production_score, na.rm = TRUE),
    complement_sd = sd(Complement_production_score, na.rm = TRUE)
  )

cat("\nAstrocyte state scores by condition:\n")
print(state_summary)

write.csv(state_summary, "output/tables/snrnaseq/astrocytes/state_scores_by_condition.csv", row.names = FALSE)

# =============================================================================
# PART 3: Statistical testing of state differences
# =============================================================================

cat("\n=== Statistical Testing ===\n")

# Compare Implant vs Control for each score
scores_to_test <- c("Astrocyte_reactive_score", "Astrocyte_homeostatic_score",
                    "A2_neuroprotective_score", "GlialScar_score", "Complement_production_score")

stat_results <- data.frame()
for(score in scores_to_test) {
  # Implant vs Control
  ctrl_vals <- astro@meta.data %>% filter(Condition == "Control") %>% pull(!!sym(score))
  impl_vals <- astro@meta.data %>% filter(Condition == "Implant") %>% pull(!!sym(score))
  stab_vals <- astro@meta.data %>% filter(Condition == "Stab") %>% pull(!!sym(score))

  # Wilcoxon tests
  impl_ctrl <- wilcox.test(impl_vals, ctrl_vals)
  stab_ctrl <- wilcox.test(stab_vals, ctrl_vals)
  impl_stab <- wilcox.test(impl_vals, stab_vals)

  stat_results <- rbind(stat_results, data.frame(
    score = score,
    comparison = c("Implant_vs_Control", "Stab_vs_Control", "Implant_vs_Stab"),
    pvalue = c(impl_ctrl$p.value, stab_ctrl$p.value, impl_stab$p.value),
    diff_mean = c(mean(impl_vals) - mean(ctrl_vals),
                  mean(stab_vals) - mean(ctrl_vals),
                  mean(impl_vals) - mean(stab_vals))
  ))
}

stat_results$fdr <- p.adjust(stat_results$pvalue, method = "BH")
stat_results$significant <- stat_results$fdr < 0.05

cat("\nSignificant differences:\n")
print(stat_results %>% filter(significant))

write.csv(stat_results, "output/tables/snrnaseq/astrocytes/state_statistics.csv", row.names = FALSE)

# =============================================================================
# PART 4: Key marker gene expression
# =============================================================================

cat("\n=== Key Marker Expression ===\n")

key_genes <- c("Gfap", "Vim", "C3", "Tnc", "Lcn2", "Serpina3n", "Aqp4", "Slc1a2", "Aldh1l1")

# Get expression data
expr_data <- FetchData(astro, vars = c(key_genes, "Condition"))

# Calculate stats per condition
marker_stats <- expr_data %>%
  pivot_longer(cols = all_of(key_genes), names_to = "gene", values_to = "expression") %>%
  group_by(Condition, gene) %>%
  summarise(
    mean_expr = mean(expression),
    pct_expressing = mean(expression > 0) * 100,
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Condition, values_from = c(mean_expr, pct_expressing))

cat("\nKey marker expression:\n")
print(marker_stats)

write.csv(marker_stats, "output/tables/snrnaseq/astrocytes/marker_expression.csv", row.names = FALSE)

# =============================================================================
# PART 5: Pseudobulk differential expression
# =============================================================================

cat("\n=== Pseudobulk Differential Expression ===\n")

# Aggregate by sample
DefaultAssay(astro) <- "RNA"

# Get sample info
sample_info <- astro@meta.data %>%
  select(orig.ident, Condition) %>%
  distinct()

# Create pseudobulk counts
pb_counts <- AggregateExpression(astro, group.by = "orig.ident", slot = "counts")$RNA

# Filter samples with enough cells
cells_per_sample <- table(astro$orig.ident)
valid_samples <- names(cells_per_sample)[cells_per_sample >= 10]
cat("Samples with >=10 astrocytes:", length(valid_samples), "\n")

if(length(valid_samples) >= 4) {
  pb_counts <- pb_counts[, valid_samples]

  # DESeq2 analysis
  library(DESeq2)

  # Get conditions for valid samples
  sample_conditions <- sample_info %>%
    filter(orig.ident %in% valid_samples) %>%
    arrange(match(orig.ident, valid_samples))

  coldata <- data.frame(
    row.names = valid_samples,
    condition = factor(sample_conditions$Condition)
  )

  # Filter low counts
  keep <- rowSums(pb_counts >= 10) >= 2
  pb_counts_filt <- pb_counts[keep, ]

  # DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = round(pb_counts_filt),
    colData = coldata,
    design = ~ condition
  )

  dds <- DESeq(dds, quiet = TRUE)

  # Implant vs Control
  res_impl <- results(dds, contrast = c("condition", "Implant", "Control"))
  res_impl_df <- as.data.frame(res_impl) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(pvalue) %>%
    mutate(comparison = "Implant_vs_Control")

  # Stab vs Control
  res_stab <- results(dds, contrast = c("condition", "Stab", "Control"))
  res_stab_df <- as.data.frame(res_stab) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(pvalue) %>%
    mutate(comparison = "Stab_vs_Control")

  # Combine
  de_results <- bind_rows(res_impl_df, res_stab_df)

  cat("\nDEGs (FDR < 0.05):\n")
  cat("Implant vs Control:", sum(res_impl_df$padj < 0.05, na.rm = TRUE), "\n")
  cat("Stab vs Control:", sum(res_stab_df$padj < 0.05, na.rm = TRUE), "\n")

  # Top genes
  cat("\nTop Implant vs Control genes:\n")
  print(head(res_impl_df %>% filter(padj < 0.05) %>% arrange(desc(abs(log2FoldChange))), 20))

  write.csv(res_impl_df, "output/tables/snrnaseq/astrocytes/de_implant_vs_control.csv", row.names = FALSE)
  write.csv(res_stab_df, "output/tables/snrnaseq/astrocytes/de_stab_vs_control.csv", row.names = FALSE)

} else {
  cat("Not enough samples for pseudobulk DE analysis\n")
}

# =============================================================================
# PART 6: Compare complement production: Astrocytes vs Microglia
# =============================================================================

cat("\n=== Astrocyte vs Microglia Complement Production ===\n")

# Get microglia too
micro <- subset(seu, celltype_l1 == "Microglia")

# Calculate C3 expression in both
c3_astro <- FetchData(astro, vars = c("C3", "Condition"))
c3_micro <- FetchData(micro, vars = c("C3", "Condition"))

c3_comparison <- data.frame(
  cell_type = c(rep("Astrocyte", nrow(c3_astro)), rep("Microglia", nrow(c3_micro))),
  C3_expr = c(c3_astro$C3, c3_micro$C3),
  Condition = c(c3_astro$Condition, c3_micro$Condition)
)

c3_summary <- c3_comparison %>%
  group_by(cell_type, Condition) %>%
  summarise(
    n_cells = n(),
    mean_C3 = mean(C3_expr),
    pct_expressing = mean(C3_expr > 0) * 100,
    .groups = "drop"
  )

cat("\nC3 expression comparison:\n")
print(c3_summary)

write.csv(c3_summary, "output/tables/snrnaseq/astrocytes/c3_astrocyte_vs_microglia.csv", row.names = FALSE)

# Also check C1q (should be mostly microglia)
c1q_astro <- FetchData(astro, vars = c("C1qa", "C1qb", "C1qc", "Condition"))
c1q_micro <- FetchData(micro, vars = c("C1qa", "C1qb", "C1qc", "Condition"))

c1q_astro$total_c1q <- c1q_astro$C1qa + c1q_astro$C1qb + c1q_astro$C1qc
c1q_micro$total_c1q <- c1q_micro$C1qa + c1q_micro$C1qb + c1q_micro$C1qc

c1q_comparison <- data.frame(
  cell_type = c(rep("Astrocyte", nrow(c1q_astro)), rep("Microglia", nrow(c1q_micro))),
  C1q_total = c(c1q_astro$total_c1q, c1q_micro$total_c1q),
  Condition = c(c1q_astro$Condition, c1q_micro$Condition)
)

c1q_summary <- c1q_comparison %>%
  group_by(cell_type, Condition) %>%
  summarise(
    n_cells = n(),
    mean_C1q = mean(C1q_total),
    pct_expressing = mean(C1q_total > 0) * 100,
    .groups = "drop"
  )

cat("\nC1q expression comparison:\n")
print(c1q_summary)

write.csv(c1q_summary, "output/tables/snrnaseq/astrocytes/c1q_astrocyte_vs_microglia.csv", row.names = FALSE)

# =============================================================================
# PART 7: Visualization
# =============================================================================

cat("\n=== Generating Figures ===\n")

# State scores violin plots
score_data <- astro@meta.data %>%
  select(Condition, Astrocyte_reactive_score, Astrocyte_homeostatic_score,
         A2_neuroprotective_score, GlialScar_score, Complement_production_score) %>%
  pivot_longer(-Condition, names_to = "score", values_to = "value") %>%
  mutate(score = gsub("_score", "", score),
         score = gsub("Astrocyte_", "", score))

p1 <- ggplot(score_data, aes(x = Condition, y = value, fill = Condition)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.size = 0.5) +
  facet_wrap(~score, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c("Control" = "#4DAF4A", "Implant" = "#E41A1C", "Stab" = "#377EB8")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Astrocyte State Scores by Condition", y = "Module Score")

ggsave("output/figures/snrnaseq/astrocytes/state_violins.png", p1, width = 14, height = 4, dpi = 300)
ggsave("output/figures/snrnaseq/astrocytes/state_violins.pdf", p1, width = 14, height = 4)

# C3 comparison plot
p2 <- ggplot(c3_comparison, aes(x = cell_type, y = C3_expr, fill = Condition)) +
  geom_violin(scale = "width", alpha = 0.7, position = position_dodge(0.8)) +
  geom_boxplot(width = 0.2, position = position_dodge(0.8), outlier.size = 0.5) +
  scale_fill_manual(values = c("Control" = "#4DAF4A", "Implant" = "#E41A1C", "Stab" = "#377EB8")) +
  theme_minimal() +
  labs(title = "C3 Expression: Astrocytes vs Microglia",
       subtitle = "Astrocytes are significant C3 producers",
       x = "Cell Type", y = "C3 Expression")

ggsave("output/figures/snrnaseq/astrocytes/c3_comparison.png", p2, width = 8, height = 5, dpi = 300)
ggsave("output/figures/snrnaseq/astrocytes/c3_comparison.pdf", p2, width = 8, height = 5)

# Key markers heatmap (remove genes with no variance)
marker_matrix <- expr_data %>%
  group_by(Condition) %>%
  summarise(across(all_of(key_genes), mean)) %>%
  tibble::column_to_rownames("Condition") %>%
  as.matrix() %>%
  t()

# Remove genes with zero variance (all zeros)
gene_var <- apply(marker_matrix, 1, var, na.rm = TRUE)
marker_matrix <- marker_matrix[gene_var > 0, , drop = FALSE]

# Scale for visualization
marker_scaled <- t(scale(t(marker_matrix)))
marker_scaled[is.na(marker_scaled)] <- 0  # Handle any remaining NAs

library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#377EB8", "white", "#E41A1C"))

ht <- Heatmap(marker_scaled,
              name = "Z-score",
              col = col_fun,
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              column_title = "Astrocyte Marker Expression by Condition")

pdf("output/figures/snrnaseq/astrocytes/marker_heatmap.pdf", width = 6, height = 6)
draw(ht)
dev.off()

png("output/figures/snrnaseq/astrocytes/marker_heatmap.png", width = 6, height = 6, units = "in", res = 300)
draw(ht)
dev.off()

# =============================================================================
# PART 8: Summary statistics for KG
# =============================================================================

cat("\n=== Summary for Knowledge Graph ===\n")

summary_stats <- list(
  total_astrocytes = ncol(astro),
  by_condition = as.list(table(astro$Condition)),
  reactive_implant_vs_ctrl = stat_results %>%
    filter(score == "Astrocyte_reactive_score", comparison == "Implant_vs_Control") %>%
    select(diff_mean, pvalue, fdr),
  complement_implant_vs_ctrl = stat_results %>%
    filter(score == "Complement_production_score", comparison == "Implant_vs_Control") %>%
    select(diff_mean, pvalue, fdr),
  c3_pct_expressing_astrocyte_implant = c3_summary %>%
    filter(cell_type == "Astrocyte", Condition == "Implant") %>%
    pull(pct_expressing),
  c3_pct_expressing_microglia_implant = c3_summary %>%
    filter(cell_type == "Microglia", Condition == "Implant") %>%
    pull(pct_expressing)
)

cat("\nKey findings:\n")
cat("1. Total astrocytes:", summary_stats$total_astrocytes, "\n")
cat("2. Reactive score Implant vs Control: diff =",
    round(summary_stats$reactive_implant_vs_ctrl$diff_mean, 4),
    ", p =", format(summary_stats$reactive_implant_vs_ctrl$pvalue, digits = 3), "\n")
cat("3. Complement production Implant vs Control: diff =",
    round(summary_stats$complement_implant_vs_ctrl$diff_mean, 4),
    ", p =", format(summary_stats$complement_implant_vs_ctrl$pvalue, digits = 3), "\n")
cat("4. C3 expressing (Implant): Astrocytes =",
    round(summary_stats$c3_pct_expressing_astrocyte_implant, 1),
    "%, Microglia =",
    round(summary_stats$c3_pct_expressing_microglia_implant, 1), "%\n")

saveRDS(summary_stats, "output/tables/snrnaseq/astrocytes/summary_stats.rds")

cat("\n=== Astrocyte analysis complete ===\n")
cat("Output tables: output/tables/snrnaseq/astrocytes/\n")
cat("Output figures: output/figures/snrnaseq/astrocytes/\n")
