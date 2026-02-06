# 06a_celltype_attribution.R
# Cell-type attribution of the 50 concordant genes using snRNA-seq
# Which cells express the implant signature?

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
source("R/config.R")

cat("=== Cell-Type Attribution of Concordant Signature ===\n")

# Create output directories
dir.create(file.path(OUT_TABLES, "snrnaseq"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_FIGURES, "manuscript"), recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. LOAD DATA
# ============================================================

cat("\n1. Loading snRNA-seq data...\n")

# Load Seurat object
seu_path <- SNRNASEQ_PATH  # defined in config.R
if (!file.exists(seu_path)) {
  stop("snRNA-seq data not found at: ", seu_path)
}

seu <- readRDS(seu_path)

# CRITICAL: Use RNA assay for expression analysis (not SCT)
DefaultAssay(seu) <- "RNA"
cat(sprintf("   Cells: %d\n", ncol(seu)))
cat(sprintf("   Genes: %d\n", nrow(seu)))
cat(sprintf("   Assay: %s\n", DefaultAssay(seu)))

# Check available metadata
cat("\n   Cell types (celltype_l2):\n")
print(table(seu$celltype_l2))

cat("\n   Conditions:\n")
print(table(seu$Condition))

# Load concordant genes
concordant <- read.csv(file.path(OUT_TABLES_COMPARISON, "validated_genes.csv"))
cat(sprintf("\n   Concordant genes: %d\n", nrow(concordant)))

# Standardize gene names (Seurat uses title case for rat)
concordant$seurat_gene <- tools::toTitleCase(tolower(concordant$gene_upper))

# Check which genes are in the Seurat object
genes_in_seu <- concordant$seurat_gene[concordant$seurat_gene %in% rownames(seu)]
cat(sprintf("   Genes found in snRNA-seq: %d/%d\n", length(genes_in_seu), nrow(concordant)))

missing_genes <- concordant$seurat_gene[!concordant$seurat_gene %in% rownames(seu)]
if (length(missing_genes) > 0) {
  cat("   Missing genes:", paste(missing_genes[1:min(10, length(missing_genes))], collapse = ", "), "\n")
}

# ============================================================
# 2. CALCULATE EXPRESSION BY CELL TYPE
# ============================================================

cat("\n2. Calculating expression by cell type...\n")

# Get normalized expression matrix for concordant genes
expr_mat <- GetAssayData(seu, slot = "data")[genes_in_seu, ]

# Get cell type annotations
cell_types <- seu$celltype_l2
conditions <- seu$Condition

# Get unique cell types
unique_celltypes <- unique(as.character(cell_types))
unique_celltypes_clean <- gsub(" ", "_", unique_celltypes)
cat("   Unique cell types: ", paste(unique_celltypes, collapse = ", "), "\n")

# Calculate mean expression per cell type
mean_expr_list <- lapply(unique_celltypes, function(ct) {
  cells <- which(cell_types == ct)
  if (length(cells) > 0) {
    rowMeans(expr_mat[, cells, drop = FALSE])
  } else {
    rep(0, length(genes_in_seu))
  }
})
mean_expr_by_celltype <- do.call(cbind, mean_expr_list)
colnames(mean_expr_by_celltype) <- unique_celltypes_clean
rownames(mean_expr_by_celltype) <- genes_in_seu

# Calculate total expression from each cell type
total_expr_list <- lapply(unique_celltypes, function(ct) {
  cells <- which(cell_types == ct)
  if (length(cells) > 0) {
    rowSums(expr_mat[, cells, drop = FALSE])
  } else {
    rep(0, length(genes_in_seu))
  }
})
total_expr_by_celltype <- do.call(cbind, total_expr_list)
colnames(total_expr_by_celltype) <- unique_celltypes_clean
rownames(total_expr_by_celltype) <- genes_in_seu

# Calculate fraction of expression from each cell type
total_expr <- rowSums(total_expr_by_celltype)
pct_expr_by_celltype <- total_expr_by_celltype / total_expr * 100
pct_expr_by_celltype[is.na(pct_expr_by_celltype)] <- 0

cat("   Cell type columns: ", paste(colnames(pct_expr_by_celltype), collapse = ", "), "\n")

# Calculate % cells expressing each gene by cell type
pct_cells_list <- lapply(unique_celltypes, function(ct) {
  cells <- which(cell_types == ct)
  if (length(cells) > 0) {
    rowMeans(expr_mat[, cells, drop = FALSE] > 0) * 100
  } else {
    rep(0, length(genes_in_seu))
  }
})
pct_cells_expressing <- do.call(cbind, pct_cells_list)
colnames(pct_cells_expressing) <- unique_celltypes_clean
rownames(pct_cells_expressing) <- genes_in_seu

cat("   Done calculating expression metrics\n")

# ============================================================
# 3. CLASSIFY GENES BY PRIMARY CELL TYPE
# ============================================================

cat("\n3. Classifying genes by primary cell type...\n")

# Define cell type categories
neuronal_types <- c("Neuron_Excitatory", "Neuron_Inhibitory")
glial_types <- c("Microglia", "Astrocyte", "Oligodendrocyte_Myelinating",
                 "Oligodendrocyte_Precursor", "Oligodendrocyte_Premyelinating")

# For each gene, calculate % neuronal vs % glial
gene_attribution <- data.frame(
  gene = genes_in_seu,
  stringsAsFactors = FALSE
)

# Add category from concordant
gene_attribution <- gene_attribution %>%
  left_join(concordant %>% dplyr::select(seurat_gene, category, lfc_impl_ctrl) %>%
              rename(gene = seurat_gene), by = "gene")

# Calculate neuronal vs glial fraction
neuronal_cols <- colnames(pct_expr_by_celltype)[colnames(pct_expr_by_celltype) %in% neuronal_types]
glial_cols <- colnames(pct_expr_by_celltype)[colnames(pct_expr_by_celltype) %in% glial_types]

gene_attribution$pct_neuronal <- rowSums(pct_expr_by_celltype[, neuronal_cols, drop = FALSE])
gene_attribution$pct_glial <- rowSums(pct_expr_by_celltype[, glial_cols, drop = FALSE])

# Find dominant cell type
gene_attribution$dominant_celltype <- colnames(pct_expr_by_celltype)[apply(pct_expr_by_celltype, 1, which.max)]
gene_attribution$dominant_pct <- apply(pct_expr_by_celltype, 1, max)

# Classify as Neuronal, Glial, or Mixed
gene_attribution$cell_class <- case_when(
  gene_attribution$pct_neuronal > 70 ~ "Neuronal",
  gene_attribution$pct_glial > 70 ~ "Glial",
  TRUE ~ "Mixed"
)

# Add specific cell type breakdown
for (ct in colnames(pct_expr_by_celltype)) {
  col_name <- paste0("pct_", gsub(" ", "_", ct))
  gene_attribution[[col_name]] <- pct_expr_by_celltype[genes_in_seu, ct]
}

# Debug: print column names
cat("   Attribution columns: ", paste(names(gene_attribution)[grepl("pct_", names(gene_attribution))], collapse = ", "), "\n")

cat("\nGene classification summary:\n")
print(table(gene_attribution$cell_class, gene_attribution$category))

# ============================================================
# 4. CONDITION-SPECIFIC EXPRESSION
# ============================================================

cat("\n4. Calculating condition-specific expression...\n")

# Calculate expression by condition AND cell type
conditions_of_interest <- c("Control", "Implant", "Stab")

condition_celltype_expr <- list()
for (cond in conditions_of_interest) {
  cond_cells <- which(conditions == cond)
  if (length(cond_cells) > 0) {
    cond_celltypes <- cell_types[cond_cells]
    cond_expr <- expr_mat[, cond_cells, drop = FALSE]

    mean_by_ct <- sapply(unique(cell_types), function(ct) {
      ct_cells <- which(cond_celltypes == ct)
      if (length(ct_cells) > 0) {
        rowMeans(cond_expr[, ct_cells, drop = FALSE])
      } else {
        rep(NA, nrow(cond_expr))
      }
    })
    rownames(mean_by_ct) <- genes_in_seu
    condition_celltype_expr[[cond]] <- mean_by_ct
  }
}

# Calculate fold change (Implant vs Control) by cell type
if ("Control" %in% names(condition_celltype_expr) && "Implant" %in% names(condition_celltype_expr)) {
  fc_by_celltype <- log2((condition_celltype_expr[["Implant"]] + 0.1) /
                          (condition_celltype_expr[["Control"]] + 0.1))

  cat("\nMean log2FC by cell type (Implant vs Control):\n")
  print(round(colMeans(fc_by_celltype, na.rm = TRUE), 2))
}

# ============================================================
# 5. SAVE RESULTS
# ============================================================

cat("\n5. Saving results...\n")

# Save gene attribution table
write.csv(gene_attribution,
          file.path(OUT_TABLES, "snrnaseq", "concordant_celltype_attribution.csv"),
          row.names = FALSE)

# Save full expression matrix by cell type
expr_summary <- as.data.frame(pct_expr_by_celltype)
expr_summary$gene <- rownames(expr_summary)
expr_summary <- expr_summary %>%
  left_join(gene_attribution %>% dplyr::select(gene, category, cell_class), by = "gene")
write.csv(expr_summary,
          file.path(OUT_TABLES, "snrnaseq", "concordant_expression_by_celltype.csv"),
          row.names = FALSE)

# ============================================================
# 6. PUBLICATION FIGURES
# ============================================================

cat("\n6. Creating publication figures...\n")

# Figure 1: Stacked bar - % expression by cell type for each gene
# Get pct columns (exclude neuronal/glial summaries)
pct_cols <- names(gene_attribution)[grepl("^pct_", names(gene_attribution))]
pct_cols <- pct_cols[!pct_cols %in% c("pct_neuronal", "pct_glial")]
cat("   Plotting columns: ", paste(pct_cols, collapse = ", "), "\n")

plot_data <- gene_attribution %>%
  dplyr::select(gene, category, all_of(pct_cols)) %>%
  pivot_longer(cols = all_of(pct_cols), names_to = "celltype", values_to = "pct") %>%
  mutate(celltype = gsub("pct_", "", celltype),
         celltype = gsub("_", " ", celltype))

# Order genes by glial fraction (most glial at top)
gene_order <- gene_attribution %>%
  arrange(desc(pct_glial)) %>%
  pull(gene)

plot_data$gene <- factor(plot_data$gene, levels = gene_order)

# Simplify cell type names
plot_data$celltype_simple <- case_when(
  grepl("Neuron", plot_data$celltype) ~ "Neurons",
  grepl("Microglia", plot_data$celltype) ~ "Microglia",
  grepl("Astrocyte", plot_data$celltype) ~ "Astrocytes",
  grepl("Oligodendrocyte", plot_data$celltype) ~ "Oligodendrocytes",
  TRUE ~ "Other"
)

# Aggregate by simplified cell type
plot_data_simple <- plot_data %>%
  group_by(gene, category, celltype_simple) %>%
  summarize(pct = sum(pct), .groups = "drop")

celltype_colors <- c(
  "Neurons" = "#4DAF4A",
  "Microglia" = "#E41A1C",
  "Astrocytes" = "#377EB8",
  "Oligodendrocytes" = "#984EA3",
  "Other" = "#BDBDBD"
)

p1 <- ggplot(plot_data_simple, aes(x = gene, y = pct, fill = celltype_simple)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = celltype_colors) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Cell-Type Attribution of Concordant Genes",
       subtitle = "% of total expression from each cell type",
       x = NULL, y = "% of Expression", fill = "Cell Type") +
  theme_publication() +
  theme(axis.text.y = element_text(size = 7),
        legend.position = "bottom")

save_figure(file.path(OUT_FIGURES, "manuscript", "fig_celltype_attribution_bar.png"),
            p1, width = 10, height = 14)

# Figure 2: Summary pie/donut - overall signature attribution
overall_attribution <- gene_attribution %>%
  summarize(
    Neurons = mean(pct_neuronal),
    Microglia = mean(pct_Microglia, na.rm = TRUE),
    Astrocytes = mean(pct_Astrocyte, na.rm = TRUE),
    Oligodendrocytes = mean(pct_Oligodendrocyte_Myelinating, na.rm = TRUE) +
                       mean(pct_Oligodendrocyte_Precursor, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "celltype", values_to = "pct")

p2 <- ggplot(overall_attribution, aes(x = "", y = pct, fill = celltype)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = celltype_colors) +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            position = position_stack(vjust = 0.5), color = "white", fontface = "bold") +
  labs(title = "Concordant Signature: Cell-Type Origin",
       subtitle = sprintf("Based on %d genes", nrow(gene_attribution)),
       fill = "Cell Type") +
  theme_void() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

save_figure(file.path(OUT_FIGURES, "manuscript", "fig_celltype_attribution_pie.png"),
            p2, width = 8, height = 6)

# Figure 3: Neuronal vs Glial scatter
p3 <- ggplot(gene_attribution, aes(x = pct_neuronal, y = pct_glial)) +
  geom_abline(intercept = 100, slope = -1, linetype = "dashed", color = COL_REF) +
  geom_point(aes(color = category, size = abs(lfc_impl_ctrl)), alpha = 0.7) +
  geom_text_repel(aes(label = gene), size = 2.5, max.overlaps = 20) +
  scale_color_manual(values = c("Implant_Up" = COL_UP, "Implant_Down" = COL_DOWN)) +
  scale_size_continuous(range = c(2, 6), name = "|log2FC|") +
  labs(title = "Concordant Genes: Neuronal vs Glial Expression",
       subtitle = "Each point = one concordant gene",
       x = "% Expression from Neurons",
       y = "% Expression from Glia",
       color = "Direction") +
  theme_publication() +
  theme(legend.position = "right")

save_figure(file.path(OUT_FIGURES, "manuscript", "fig_celltype_neuronal_vs_glial.png"),
            p3, width = 10, height = 8)

# Figure 4: Heatmap of expression by cell type
heatmap_mat <- as.matrix(pct_expr_by_celltype)

# Scale for visualization
heatmap_scaled <- t(scale(t(heatmap_mat)))
heatmap_scaled[heatmap_scaled > 2] <- 2
heatmap_scaled[heatmap_scaled < -2] <- -2

heatmap_df <- as.data.frame(heatmap_scaled) %>%
  mutate(gene = rownames(heatmap_scaled)) %>%
  left_join(gene_attribution %>% dplyr::select(gene, category, cell_class), by = "gene") %>%
  pivot_longer(cols = -c(gene, category, cell_class), names_to = "celltype", values_to = "zscore")

# Order by cell class then gene
gene_order_hm <- gene_attribution %>%
  arrange(cell_class, desc(pct_glial)) %>%
  pull(gene)
heatmap_df$gene <- factor(heatmap_df$gene, levels = gene_order_hm)

p4 <- ggplot(heatmap_df, aes(x = celltype, y = gene, fill = zscore)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Z-score",
                       limits = c(-2, 2)) +
  facet_grid(cell_class ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Concordant Gene Expression by Cell Type",
       subtitle = "Row-scaled % expression",
       x = NULL, y = NULL) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 6),
        strip.text.y = element_text(angle = 0))

save_figure(file.path(OUT_FIGURES, "manuscript", "fig_celltype_heatmap.png"),
            p4, width = 10, height = 14)

# Figure 5: Key finding - UP genes are glial, DOWN genes are neuronal
summary_by_direction <- gene_attribution %>%
  group_by(category) %>%
  summarize(
    n = n(),
    mean_neuronal = mean(pct_neuronal),
    mean_glial = mean(pct_glial),
    se_neuronal = sd(pct_neuronal) / sqrt(n()),
    se_glial = sd(pct_glial) / sqrt(n()),
    .groups = "drop"
  )

summary_long <- summary_by_direction %>%
  pivot_longer(cols = c(mean_neuronal, mean_glial),
               names_to = "source", values_to = "pct") %>%
  mutate(source = gsub("mean_", "", source),
         source = tools::toTitleCase(source),
         se = ifelse(source == "Neuronal", se_neuronal, se_glial))

p5 <- ggplot(summary_long, aes(x = category, y = pct, fill = source)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pct - se, ymax = pct + se),
                position = position_dodge(0.7), width = 0.2) +
  scale_fill_manual(values = c("Neuronal" = "#4DAF4A", "Glial" = "#E41A1C")) +
  labs(title = "Cell-Type Origin by Gene Direction",
       subtitle = "Upregulated genes = Glial; Downregulated genes = Neuronal",
       x = NULL, y = "% of Expression", fill = "Source") +
  theme_publication() +
  theme(legend.position = "bottom")

save_figure(file.path(OUT_FIGURES, "manuscript", "fig_celltype_by_direction.png"),
            p5, width = 7, height = 6)

# ============================================================
# 7. SUMMARY STATISTICS
# ============================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("CELL-TYPE ATTRIBUTION SUMMARY\n")
cat(strrep("=", 60), "\n")

cat("\n[1] OVERALL ATTRIBUTION:\n")
cat(sprintf("    Neuronal: %.1f%%\n", mean(gene_attribution$pct_neuronal)))
cat(sprintf("    Glial: %.1f%%\n", mean(gene_attribution$pct_glial)))
cat(sprintf("    - Microglia: %.1f%%\n", mean(gene_attribution$pct_Microglia, na.rm = TRUE)))
cat(sprintf("    - Astrocytes: %.1f%%\n", mean(gene_attribution$pct_Astrocyte, na.rm = TRUE)))

cat("\n[2] BY DIRECTION:\n")
up_genes <- gene_attribution %>% filter(category == "Implant_Up")
down_genes <- gene_attribution %>% filter(category == "Implant_Down")

cat(sprintf("    UP genes (n=%d):\n", nrow(up_genes)))
cat(sprintf("      Neuronal: %.1f%%\n", mean(up_genes$pct_neuronal)))
cat(sprintf("      Glial: %.1f%%\n", mean(up_genes$pct_glial)))

cat(sprintf("    DOWN genes (n=%d):\n", nrow(down_genes)))
cat(sprintf("      Neuronal: %.1f%%\n", mean(down_genes$pct_neuronal)))
cat(sprintf("      Glial: %.1f%%\n", mean(down_genes$pct_glial)))

cat("\n[3] KEY GENES:\n")
key_genes_attr <- gene_attribution %>%
  filter(gene %in% c("Npas4", "C1qa", "C1qb", "C1qc", "C3", "Gfap", "Spp1", "Tyrobp"))

for (i in 1:nrow(key_genes_attr)) {
  g <- key_genes_attr[i, ]
  cat(sprintf("    %s: %.0f%% %s (dominant: %s)\n",
              g$gene, g$dominant_pct, g$cell_class, g$dominant_celltype))
}

cat("\n[4] STATISTICAL TEST:\n")
# Test if UP genes are more glial than DOWN genes
if (nrow(up_genes) > 2 && nrow(down_genes) >= 1) {
  wilcox_test <- wilcox.test(up_genes$pct_glial, down_genes$pct_glial, alternative = "greater")
  cat(sprintf("    UP vs DOWN glial %%: Wilcoxon p = %.2e\n", wilcox_test$p.value))
}

cat("\nSaved results to:\n")
cat(sprintf("  %s/snrnaseq/concordant_celltype_attribution.csv\n", OUT_TABLES))
cat(sprintf("  %s/manuscript/fig_celltype_*.png (5 figures)\n", OUT_FIGURES))
