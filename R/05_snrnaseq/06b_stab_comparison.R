# 06b_stab_comparison.R
# Compare concordant genes in stab wound vs implant
# Uses GSE226211 stab wound scRNA-seq to validate implant-specificity
# Also analyzes NPAS4 distance gradient from silicon timecourse

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
source("R/config.R")

cat("=== Stab Wound Comparison & NPAS4 Gradient Analysis ===\n")

# Create output directories
dir.create(file.path(OUT_TABLES, "stab_comparison"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_FIGURES, "manuscript"), recursive = TRUE, showWarnings = FALSE)

# ============================================================
# PART 1: NPAS4 DISTANCE GRADIENT (Silicon Timecourse)
# ============================================================

cat("\n=== PART 1: NPAS4 Distance Gradient ===\n")

# Load silicon timecourse DE results - extract NPAS4
tc_files <- list.files(file.path(OUT_TABLES, "timecourse"),
                       pattern = "^de_.*\\.csv$", full.names = TRUE)

npas4_list <- list()
for (f in tc_files) {
  df <- read.csv(f)
  # Find NPAS4 (case insensitive)
  npas4_idx <- which(toupper(df$gene) == "NPAS4")
  if (length(npas4_idx) > 0) {
    npas4 <- df[npas4_idx[1], ]
    comp <- gsub("^de_|\\.csv$", "", basename(f))
    parts <- strsplit(comp, "_")[[1]]
    npas4_list[[length(npas4_list) + 1]] <- data.frame(
      timepoint = parts[1],
      distance = parts[2],
      log2FC = as.numeric(npas4$log2FoldChange),
      padj = as.numeric(npas4$padj),
      significant = as.numeric(npas4$padj) < 0.05,
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %s: NPAS4 log2FC = %.2f, padj = %.2e\n", comp, npas4$log2FoldChange, npas4$padj))
  }
}
npas4_data <- do.call(rbind, npas4_list)

cat("\nNPAS4 by Distance and Time:\n")
print(npas4_data %>% arrange(timepoint, distance))

# Save NPAS4 gradient data
write.csv(npas4_data, file.path(OUT_TABLES, "stab_comparison", "npas4_distance_gradient.csv"), row.names = FALSE)

# Figure: NPAS4 gradient
npas4_data$timepoint <- factor(npas4_data$timepoint, levels = c("24h", "1wk", "6wk"))
npas4_data$distance <- factor(npas4_data$distance, levels = c("100um", "500um"))

p_npas4 <- ggplot(npas4_data, aes(x = timepoint, y = log2FC, fill = distance)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = ifelse(significant, "*", ""),
                y = log2FC - sign(log2FC) * 2),
            position = position_dodge(0.7), size = 6) +
  scale_fill_manual(values = c("100um" = "#E41A1C", "500um" = "#377EB8")) +
  labs(title = "NPAS4 Silencing: Distance Gradient",
       subtitle = "100μm = kill zone; 500μm = protected; * = FDR < 0.05",
       x = "Time Post-Implant", y = "NPAS4 log2FC (vs Naive)",
       fill = "Distance") +
  theme_publication() +
  theme(legend.position = "bottom")

save_figure(file.path(OUT_FIGURES, "manuscript", "fig_npas4_distance_gradient.png"),
            p_npas4, width = 8, height = 6)

# ============================================================
# PART 2: LOAD SNRNA-SEQ REFERENCE (has Control, Implant, Stab)
# ============================================================

cat("\n=== PART 2: Loading snRNA-seq Reference ===\n")

# Use the main snRNA-seq reference which already has Control, Implant, and Stab conditions
seu_path <- SNRNASEQ_PATH  # defined in config.R
seu <- readRDS(seu_path)

cat(sprintf("Loaded: %d cells, %d genes\n", ncol(seu), nrow(seu)))
cat("Conditions:\n")
print(table(seu$Condition))

# ============================================================
# PART 3: PSEUDOBULK COMPARISON
# ============================================================

cat("\n=== PART 3: Pseudobulk DE Analysis ===\n")

# Load concordant genes
concordant <- read.csv(file.path(OUT_TABLES_COMPARISON, "validated_genes.csv"))
concordant_genes <- unique(toupper(concordant$gene_upper))
cat(sprintf("Concordant genes to test: %d\n", length(concordant_genes)))

# Standardize gene names in Seurat (title case for rat)
concordant_seurat <- tools::toTitleCase(tolower(concordant_genes))
genes_found <- concordant_seurat[concordant_seurat %in% rownames(seu)]
cat(sprintf("Genes found in scRNA-seq: %d\n", length(genes_found)))

# Create pseudobulk by sample
# Get raw counts from RNA assay (not SCT which doesn't support JoinLayers)
cat("Available assays:", paste(names(seu@assays), collapse = ", "), "\n")
cat("Default assay:", DefaultAssay(seu), "\n")

# Switch to RNA assay for raw counts
if ("RNA" %in% names(seu@assays)) {
  DefaultAssay(seu) <- "RNA"
  cat("Switched to RNA assay for pseudobulk\n")
}

# Handle Seurat v5 multiple layers
if (packageVersion("SeuratObject") >= "5.0.0") {
  # Check if layers need joining
  rna_assay <- seu[["RNA"]]
  if (inherits(rna_assay, "Assay5")) {
    seu[["RNA"]] <- JoinLayers(rna_assay)
    cat("Joined RNA layers for Seurat v5\n")
  }
}
counts <- GetAssayData(seu, layer = "counts")

# Aggregate by sample (orig.ident)
samples <- unique(seu$orig.ident)
pseudobulk <- sapply(samples, function(s) {
  cells <- which(seu$orig.ident == s)
  if (length(cells) > 10) {
    rowSums(counts[, cells])
  } else {
    rep(NA, nrow(counts))
  }
})
rownames(pseudobulk) <- rownames(counts)
pseudobulk <- pseudobulk[, !is.na(pseudobulk[1, ])]

# Get condition for each sample (handle both Condition and condition)
cond_col <- if ("Condition" %in% colnames(seu@meta.data)) "Condition" else "condition"
sample_condition <- sapply(colnames(pseudobulk), function(s) {
  unique(seu@meta.data[[cond_col]][seu$orig.ident == s])[1]
})

cat("\nPseudobulk samples by condition:\n")
print(table(sample_condition))

# ============================================================
# PART 4: STAB vs CONTROL DE
# ============================================================

cat("\n=== PART 4: Stab vs Control DE ===\n")

# Filter to Stab and Control only
stab_ctrl_idx <- sample_condition %in% c("Stab", "Control")
if (sum(stab_ctrl_idx) < 4) {
  cat("Not enough Stab/Control samples for DE. Using available data...\n")
} else {
  pb_stab <- pseudobulk[, stab_ctrl_idx]
  cond_stab <- factor(sample_condition[stab_ctrl_idx], levels = c("Control", "Stab"))

  # DESeq2
  dds_stab <- DESeqDataSetFromMatrix(
    countData = round(pb_stab),
    colData = data.frame(condition = cond_stab),
    design = ~ condition
  )
  dds_stab <- dds_stab[rowSums(counts(dds_stab)) > 10, ]
  dds_stab <- DESeq(dds_stab)
  res_stab <- results(dds_stab, contrast = c("condition", "Stab", "Control"))
  res_stab_df <- as.data.frame(res_stab)
  res_stab_df$gene <- rownames(res_stab_df)

  cat(sprintf("Stab vs Control: %d genes tested\n", nrow(res_stab_df)))
  cat(sprintf("Significant (FDR < 0.05): %d\n", sum(res_stab_df$padj < 0.05, na.rm = TRUE)))

  write.csv(res_stab_df, file.path(OUT_TABLES, "stab_comparison", "stab_vs_control_de.csv"), row.names = FALSE)
}

# ============================================================
# PART 5: IMPLANT vs CONTROL DE (for comparison)
# ============================================================

cat("\n=== PART 5: Implant vs Control DE ===\n")

impl_ctrl_idx <- sample_condition %in% c("Implant", "Control")
if (sum(impl_ctrl_idx) >= 4) {
  pb_impl <- pseudobulk[, impl_ctrl_idx]
  cond_impl <- factor(sample_condition[impl_ctrl_idx], levels = c("Control", "Implant"))

  dds_impl <- DESeqDataSetFromMatrix(
    countData = round(pb_impl),
    colData = data.frame(condition = cond_impl),
    design = ~ condition
  )
  dds_impl <- dds_impl[rowSums(counts(dds_impl)) > 10, ]
  dds_impl <- DESeq(dds_impl)
  res_impl <- results(dds_impl, contrast = c("condition", "Implant", "Control"))
  res_impl_df <- as.data.frame(res_impl)
  res_impl_df$gene <- rownames(res_impl_df)

  cat(sprintf("Implant vs Control: %d genes tested\n", nrow(res_impl_df)))
  cat(sprintf("Significant (FDR < 0.05): %d\n", sum(res_impl_df$padj < 0.05, na.rm = TRUE)))

  write.csv(res_impl_df, file.path(OUT_TABLES, "stab_comparison", "implant_vs_control_de.csv"), row.names = FALSE)
}

# ============================================================
# PART 6: CLASSIFY CONCORDANT GENES
# ============================================================

cat("\n=== PART 6: Classifying Concordant Genes ===\n")

if (exists("res_stab_df") && exists("res_impl_df")) {
  # Match concordant genes
  res_stab_df$gene_upper <- toupper(res_stab_df$gene)
  res_impl_df$gene_upper <- toupper(res_impl_df$gene)

  cat("Concordant columns:", paste(colnames(concordant), collapse = ", "), "\n")
  cat("res_impl_df columns:", paste(colnames(res_impl_df), collapse = ", "), "\n")

  # Prepare join tables with proper renaming
  impl_for_join <- res_impl_df[, c("gene_upper", "log2FoldChange", "padj")]
  colnames(impl_for_join) <- c("gene_upper", "sc_impl_lfc", "sc_impl_padj")

  stab_for_join <- res_stab_df[, c("gene_upper", "log2FoldChange", "padj")]
  colnames(stab_for_join) <- c("gene_upper", "sc_stab_lfc", "sc_stab_padj")

  # Check if lfc_impl_ctrl exists
  conc_cols <- c("gene_upper", "category")
  if ("lfc_impl_ctrl" %in% colnames(concordant)) {
    conc_cols <- c(conc_cols, "lfc_impl_ctrl")
  }

  classification <- concordant[, conc_cols] %>%
    left_join(impl_for_join, by = "gene_upper") %>%
    left_join(stab_for_join, by = "gene_upper")

  # Classify
  classification <- classification %>%
    mutate(
      sig_impl = sc_impl_padj < 0.05 & !is.na(sc_impl_padj),
      sig_stab = sc_stab_padj < 0.05 & !is.na(sc_stab_padj),
      impl_dir = sign(sc_impl_lfc),
      stab_dir = sign(sc_stab_lfc),

      specificity = case_when(
        sig_impl & !sig_stab ~ "Implant-Specific",
        sig_impl & sig_stab & impl_dir == stab_dir ~ "Shared Injury",
        sig_impl & sig_stab & impl_dir != stab_dir ~ "Opposite",
        !sig_impl & sig_stab ~ "Stab-Only",
        TRUE ~ "Neither Significant"
      )
    )

  cat("\nGene Classification (scRNA-seq pseudobulk):\n")
  print(table(classification$specificity, classification$category))

  # NPAS4 specifically
  npas4_class <- classification[classification$gene_upper == "NPAS4", ]
  if (nrow(npas4_class) > 0) {
    cat("\nNPAS4 in scRNA-seq pseudobulk:\n")
    cat(sprintf("  Implant log2FC: %.2f (p = %.2e)\n", npas4_class$sc_impl_lfc, npas4_class$sc_impl_padj))
    cat(sprintf("  Stab log2FC: %.2f (p = %.2e)\n", npas4_class$sc_stab_lfc, npas4_class$sc_stab_padj))
    cat(sprintf("  Classification: %s\n", npas4_class$specificity))
  }

  write.csv(classification, file.path(OUT_TABLES, "stab_comparison", "concordant_gene_classification.csv"), row.names = FALSE)

  # ============================================================
  # PART 7: FIGURES
  # ============================================================

  cat("\n=== PART 7: Creating Figures ===\n")

  # Figure: Implant vs Stab log2FC for concordant genes
  plot_data <- classification %>%
    filter(!is.na(sc_impl_lfc) & !is.na(sc_stab_lfc))

  p_scatter <- ggplot(plot_data, aes(x = sc_stab_lfc, y = sc_impl_lfc)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = COL_REF) +
    geom_hline(yintercept = 0, color = "grey50") +
    geom_vline(xintercept = 0, color = "grey50") +
    geom_point(aes(color = specificity), size = 3, alpha = 0.7) +
    geom_text_repel(aes(label = gene_upper), size = 2.5, max.overlaps = 15) +
    scale_color_manual(values = c(
      "Implant-Specific" = "#E41A1C",
      "Shared Injury" = "#FF7F00",
      "Stab-Only" = "#377EB8",
      "Neither Significant" = "#BDBDBD",
      "Opposite" = "#984EA3"
    )) +
    labs(title = "Concordant Genes: Implant vs Stab Response",
         subtitle = "scRNA-seq pseudobulk (snRNA-seq reference)",
         x = "Stab vs Control log2FC",
         y = "Implant vs Control log2FC",
         color = "Classification") +
    theme_publication() +
    theme(legend.position = "right")

  save_figure(file.path(OUT_FIGURES, "manuscript", "fig_implant_vs_stab_scatter.png"),
              p_scatter, width = 10, height = 8)

  # Figure: Classification summary bar
  class_summary <- classification %>%
    group_by(category, specificity) %>%
    summarize(n = n(), .groups = "drop")

  p_bar <- ggplot(class_summary, aes(x = category, y = n, fill = specificity)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c(
      "Implant-Specific" = "#E41A1C",
      "Shared Injury" = "#FF7F00",
      "Stab-Only" = "#377EB8",
      "Neither Significant" = "#BDBDBD",
      "Opposite" = "#984EA3"
    )) +
    labs(title = "Concordant Gene Classification",
         subtitle = "Implant-specific vs Shared injury (scRNA-seq validation)",
         x = NULL, y = "Number of Genes",
         fill = "Classification") +
    theme_publication() +
    theme(legend.position = "right")

  save_figure(file.path(OUT_FIGURES, "manuscript", "fig_gene_classification_bar.png"),
              p_bar, width = 9, height = 6)
}

# ============================================================
# SUMMARY
# ============================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("STAB COMPARISON SUMMARY\n")
cat(strrep("=", 60), "\n")

cat("\n[1] NPAS4 DISTANCE GRADIENT (Silicon):\n")
for (i in 1:nrow(npas4_data)) {
  row <- npas4_data[i, ]
  sig <- ifelse(row$significant, "***", "NS")
  cat(sprintf("    %s %s: log2FC = %.1f (%s)\n", row$timepoint, row$distance, row$log2FC, sig))
}

if (exists("classification")) {
  cat("\n[2] CONCORDANT GENE CLASSIFICATION:\n")
  for (spec in unique(classification$specificity)) {
    n <- sum(classification$specificity == spec)
    cat(sprintf("    %s: %d genes\n", spec, n))
  }

  impl_spec <- sum(classification$specificity == "Implant-Specific")
  shared <- sum(classification$specificity == "Shared Injury")
  cat(sprintf("\n[3] IMPLANT-SPECIFICITY RATE: %d/%d (%.1f%%) truly implant-specific\n",
              impl_spec, impl_spec + shared, 100 * impl_spec / (impl_spec + shared)))
}

cat("\nSaved results to:\n")
cat(sprintf("  %s/stab_comparison/\n", OUT_TABLES))
cat(sprintf("  %s/manuscript/fig_npas4_distance_gradient.png\n", OUT_FIGURES))
