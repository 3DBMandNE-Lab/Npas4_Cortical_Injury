# 06h_spp1_microglia_analysis.R
# Deep characterization of SPP1+ microglia biology
# Compare SPP1-high vs SPP1-low microglia, run pathway enrichment

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DESeq2)
library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
source("R/config.R")

cat("=== SPP1+ Microglia Characterization ===\n\n")

# Create output directory
dir.create(file.path(OUT_TABLES_SNRNASEQ, "spp1_microglia"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_FIGURES_SNRNASEQ, "spp1_microglia"), recursive = TRUE, showWarnings = FALSE)

# Load snRNA-seq data
SNRNASEQ_SOURCE <- SNRNASEQ_PATH  # defined in config.R
cat("Loading snRNA-seq data...\n")
seu <- readRDS(SNRNASEQ_SOURCE)
cat(sprintf("Loaded: %d cells, %d genes\n", ncol(seu), nrow(seu)))

# Use RNA assay
DefaultAssay(seu) <- "RNA"

# Find microglia
celltype_col <- if ("celltype_l2" %in% colnames(seu@meta.data)) "celltype_l2" else "celltype_l1"
cat(sprintf("Using cell type column: %s\n", celltype_col))

# Identify microglia cells
microglia_patterns <- c("Microglia", "microglia", "MG")
all_celltypes <- unique(seu@meta.data[[celltype_col]])
microglia_types <- all_celltypes[grepl(paste(microglia_patterns, collapse = "|"), all_celltypes, ignore.case = TRUE)]

cat(sprintf("Microglia cell types found: %s\n", paste(microglia_types, collapse = ", ")))

# Subset to microglia
seu_mg <- subset(seu, cells = colnames(seu)[seu@meta.data[[celltype_col]] %in% microglia_types])
DefaultAssay(seu_mg) <- "RNA"  # CRITICAL: Ensure RNA assay after subset
cat(sprintf("Microglia subset: %d cells\n", ncol(seu_mg)))

# Check SPP1 expression
if (!"Spp1" %in% rownames(seu_mg)) {
  stop("Spp1 not found in gene names")
}

# Get SPP1 expression - CRITICAL: Use RNA assay explicitly
spp1_expr <- GetAssayData(seu_mg, assay = "RNA", layer = "data")["Spp1", ]
cat(sprintf("SPP1 expression range: %.2f - %.2f\n", min(spp1_expr), max(spp1_expr)))
cat(sprintf("SPP1+ cells (>0): %d (%.1f%%)\n", sum(spp1_expr > 0), 100 * mean(spp1_expr > 0)))

# Define SPP1+ vs SPP1- (binary since SPP1 is sparsely expressed)
# SPP1+ = any expression, SPP1- = zero expression
seu_mg$spp1_group <- ifelse(spp1_expr > 0, "SPP1_Pos", "SPP1_Neg")

# Also create a continuous score for later analysis
seu_mg$spp1_expr <- spp1_expr

cat("\nSPP1 group distribution:\n")
print(table(seu_mg$spp1_group))

cat("\nSPP1 groups by condition:\n")
print(table(seu_mg$spp1_group, seu_mg$Condition))

# ============================================================
# PART 1: SPP1+ enrichment by condition
# ============================================================

cat("\n=== PART 1: SPP1+ Enrichment by Condition ===\n")

spp1_by_cond <- seu_mg@meta.data %>%
  group_by(Condition) %>%
  summarize(
    n_cells = n(),
    n_spp1_high = sum(spp1_group == "SPP1_Pos"),
    pct_spp1_high = 100 * n_spp1_high / n_cells,
    mean_spp1 = mean(spp1_expr[match(row.names(cur_data()), names(spp1_expr))]),
    .groups = "drop"
  )

# Recalculate mean SPP1 properly
for (cond in unique(seu_mg$Condition)) {
  cells <- colnames(seu_mg)[seu_mg$Condition == cond]
  spp1_by_cond$mean_spp1[spp1_by_cond$Condition == cond] <- mean(spp1_expr[cells])
}

cat("\nSPP1 expression by condition:\n")
print(spp1_by_cond)

write.csv(spp1_by_cond, file.path(OUT_TABLES_SNRNASEQ, "spp1_microglia", "spp1_by_condition.csv"), row.names = FALSE)

# ============================================================
# PART 2: Differential Expression SPP1-High vs SPP1-Low
# ============================================================

cat("\n=== PART 2: DE Analysis SPP1-High vs SPP1-Low ===\n")

# Filter to High and Low only for cleaner comparison
seu_mg_hl <- subset(seu_mg, cells = colnames(seu_mg)[seu_mg$spp1_group %in% c("SPP1_Pos", "SPP1_Neg")])
cat(sprintf("Cells for DE: %d SPP1_Pos, %d SPP1_Neg\n",
            sum(seu_mg_hl$spp1_group == "SPP1_Pos"),
            sum(seu_mg_hl$spp1_group == "SPP1_Neg")))

# Run Seurat FindMarkers (faster than pseudobulk for this)
Idents(seu_mg_hl) <- "spp1_group"
markers <- FindMarkers(seu_mg_hl, ident.1 = "SPP1_Pos", ident.2 = "SPP1_Neg",
                       min.pct = 0.1, logfc.threshold = 0.25, test.use = "wilcox")

markers$gene <- rownames(markers)
markers <- markers %>% arrange(p_val_adj, desc(abs(avg_log2FC)))

cat(sprintf("\nDE results: %d genes tested\n", nrow(markers)))
cat(sprintf("Significant (FDR < 0.05): %d\n", sum(markers$p_val_adj < 0.05)))
cat(sprintf("  Upregulated in SPP1_Pos: %d\n", sum(markers$p_val_adj < 0.05 & markers$avg_log2FC > 0)))
cat(sprintf("  Downregulated in SPP1_Pos: %d\n", sum(markers$p_val_adj < 0.05 & markers$avg_log2FC < 0)))

write.csv(markers, file.path(OUT_TABLES_SNRNASEQ, "spp1_microglia", "de_spp1high_vs_low.csv"), row.names = FALSE)

# Top genes
cat("\nTop 20 genes upregulated in SPP1_Pos microglia:\n")
top_up <- markers %>% filter(p_val_adj < 0.05, avg_log2FC > 0) %>% head(20)
print(top_up[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")])

cat("\nTop 20 genes downregulated in SPP1_Pos microglia:\n")
top_down <- markers %>% filter(p_val_adj < 0.05, avg_log2FC < 0) %>% head(20)
print(top_down[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")])

# ============================================================
# PART 3: GO Enrichment Analysis
# ============================================================

cat("\n=== PART 3: GO Enrichment Analysis ===\n")

# Get significant genes
sig_up <- markers %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>% pull(gene)
sig_down <- markers %>% filter(p_val_adj < 0.05, avg_log2FC < -0.5) %>% pull(gene)

cat(sprintf("Genes for enrichment: %d up, %d down (|log2FC| > 0.5)\n", length(sig_up), length(sig_down)))

# Convert to Entrez IDs
convert_to_entrez <- function(genes) {
  # Try title case (rat genes)
  genes_title <- tools::toTitleCase(tolower(genes))
  mapped <- bitr(genes_title, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)
  return(mapped$ENTREZID)
}

# GO enrichment for upregulated genes
if (length(sig_up) >= 10) {
  entrez_up <- convert_to_entrez(sig_up)
  cat(sprintf("Mapped %d/%d upregulated genes to Entrez IDs\n", length(entrez_up), length(sig_up)))

  if (length(entrez_up) >= 5) {
    ego_up_bp <- enrichGO(gene = entrez_up,
                          OrgDb = org.Rn.eg.db,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.1,
                          readable = TRUE)

    if (!is.null(ego_up_bp) && nrow(ego_up_bp) > 0) {
      cat(sprintf("\nGO BP enriched (SPP1_Pos UP): %d terms\n", nrow(ego_up_bp)))

      ego_up_df <- as.data.frame(ego_up_bp)
      write.csv(ego_up_df, file.path(OUT_TABLES_SNRNASEQ, "spp1_microglia", "go_bp_spp1high_up.csv"), row.names = FALSE)

      cat("\nTop GO BP terms (SPP1_Pos upregulated):\n")
      print(head(ego_up_df[, c("Description", "GeneRatio", "p.adjust", "geneID")], 15))

      # Plot
      if (nrow(ego_up_bp) >= 3) {
        p_go_up <- dotplot(ego_up_bp, showCategory = 20, title = "GO BP: SPP1-High Microglia (Upregulated)")
        ggsave(file.path(OUT_FIGURES_SNRNASEQ, "spp1_microglia", "go_bp_spp1high_up.png"),
               p_go_up, width = 10, height = 8)
      }
    }

    # Cellular Component
    ego_up_cc <- enrichGO(gene = entrez_up,
                          OrgDb = org.Rn.eg.db,
                          ont = "CC",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.1,
                          readable = TRUE)

    if (!is.null(ego_up_cc) && nrow(ego_up_cc) > 0) {
      ego_up_cc_df <- as.data.frame(ego_up_cc)
      write.csv(ego_up_cc_df, file.path(OUT_TABLES_SNRNASEQ, "spp1_microglia", "go_cc_spp1high_up.csv"), row.names = FALSE)
      cat(sprintf("\nGO CC enriched (SPP1_Pos UP): %d terms\n", nrow(ego_up_cc)))
    }
  }
}

# GO enrichment for downregulated genes
if (length(sig_down) >= 10) {
  entrez_down <- convert_to_entrez(sig_down)
  cat(sprintf("\nMapped %d/%d downregulated genes to Entrez IDs\n", length(entrez_down), length(sig_down)))

  if (length(entrez_down) >= 5) {
    ego_down_bp <- enrichGO(gene = entrez_down,
                            OrgDb = org.Rn.eg.db,
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.1,
                            readable = TRUE)

    if (!is.null(ego_down_bp) && nrow(ego_down_bp) > 0) {
      cat(sprintf("\nGO BP enriched (SPP1_Pos DOWN): %d terms\n", nrow(ego_down_bp)))

      ego_down_df <- as.data.frame(ego_down_bp)
      write.csv(ego_down_df, file.path(OUT_TABLES_SNRNASEQ, "spp1_microglia", "go_bp_spp1high_down.csv"), row.names = FALSE)

      cat("\nTop GO BP terms (SPP1_Pos downregulated):\n")
      print(head(ego_down_df[, c("Description", "GeneRatio", "p.adjust", "geneID")], 15))

      if (nrow(ego_down_bp) >= 3) {
        p_go_down <- dotplot(ego_down_bp, showCategory = 20, title = "GO BP: SPP1-High Microglia (Downregulated)")
        ggsave(file.path(OUT_FIGURES_SNRNASEQ, "spp1_microglia", "go_bp_spp1high_down.png"),
               p_go_down, width = 10, height = 8)
      }
    }
  }
}

# ============================================================
# PART 4: KEGG Pathway Enrichment
# ============================================================

cat("\n=== PART 4: KEGG Pathway Enrichment ===\n")

if (length(sig_up) >= 10 && exists("entrez_up") && length(entrez_up) >= 5) {
  kegg_up <- enrichKEGG(gene = entrez_up,
                        organism = "rno",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.1)

  if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
    kegg_up_df <- as.data.frame(kegg_up)
    write.csv(kegg_up_df, file.path(OUT_TABLES_SNRNASEQ, "spp1_microglia", "kegg_spp1high_up.csv"), row.names = FALSE)
    cat(sprintf("\nKEGG enriched (SPP1_Pos UP): %d pathways\n", nrow(kegg_up)))
    cat("\nTop KEGG pathways:\n")
    print(head(kegg_up_df[, c("Description", "GeneRatio", "p.adjust")], 10))
  }
}

# ============================================================
# PART 5: Signature Analysis - Known Microglial States
# ============================================================

cat("\n=== PART 5: Known Microglial State Signatures ===\n")

# Define comprehensive microglial signatures
microglial_signatures <- list(
  # Core states
  Homeostatic = c("P2ry12", "Tmem119", "Cx3cr1", "Siglech", "Hexb", "Cst3", "Sall1"),
  DAM_Stage1 = c("Tyrobp", "Ctsd", "Lyz2", "Ctsb", "Apoe"),
  DAM_Stage2 = c("Trem2", "Axl", "Cst7", "Lpl", "Cd9", "Csf1", "Itgax", "Clec7a", "Lilrb4"),

  # Functional states
  Phagocytic = c("Cd68", "Mertk", "Axl", "Tyro3", "Megf10", "Gulp1"),
  Antigen_Presenting = c("Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1", "Ciita"),
  Inflammatory_M1 = c("Il1b", "Tnf", "Il6", "Nos2", "Cd86", "Ccl2", "Cxcl10"),
  Anti_inflammatory_M2 = c("Arg1", "Mrc1", "Cd163", "Il10", "Tgfb1", "Chil3"),

  # Proliferative/activated
  Proliferative = c("Mki67", "Top2a", "Pcna", "Cdk1", "Ccnb1"),
  IFN_Response = c("Ifit1", "Ifit2", "Ifit3", "Isg15", "Mx1", "Stat1", "Irf7"),

  # Lipid metabolism (relevant to SPP1)
  Lipid_Metabolism = c("Apoe", "Lpl", "Fabp5", "Abca1", "Abcg1", "Plin2"),

  # Complement
  Complement_Production = c("C1qa", "C1qb", "C1qc", "C3", "C4b"),

  # Synapse pruning
  Synapse_Pruning = c("C1qa", "C1qb", "C1qc", "C3", "Trem2", "Tyrobp"),

  # SPP1-associated (from literature)
  SPP1_Module = c("Spp1", "Gpnmb", "Lgals3", "Fabp5", "Igf1", "Cd63", "Lpl")
)

# Check gene availability and score
available_genes <- rownames(seu_mg)
cat("\nSignature gene availability:\n")
sig_availability <- data.frame()

for (sig_name in names(microglial_signatures)) {
  genes <- microglial_signatures[[sig_name]]
  # Convert to title case for rat
  genes_title <- tools::toTitleCase(tolower(genes))
  present <- intersect(genes_title, available_genes)

  cat(sprintf("  %s: %d/%d genes\n", sig_name, length(present), length(genes)))
  sig_availability <- rbind(sig_availability, data.frame(
    signature = sig_name,
    n_genes = length(genes),
    n_present = length(present),
    genes_present = paste(present, collapse = ", ")
  ))

  # Score if enough genes
  if (length(present) >= 3) {
    seu_mg <- AddModuleScore(seu_mg, features = list(present), name = sig_name, seed = 42)
    colnames(seu_mg@meta.data)[ncol(seu_mg@meta.data)] <- sig_name
  }
}

write.csv(sig_availability, file.path(OUT_TABLES_SNRNASEQ, "spp1_microglia", "signature_availability.csv"), row.names = FALSE)

# Compare signature scores between SPP1 groups
sig_cols <- intersect(names(microglial_signatures), colnames(seu_mg@meta.data))
cat(sprintf("\n%d signatures scored\n", length(sig_cols)))

if (length(sig_cols) > 0) {
  sig_comparison <- seu_mg@meta.data %>%
    filter(spp1_group %in% c("SPP1_Pos", "SPP1_Neg")) %>%
    group_by(spp1_group) %>%
    summarize(across(all_of(sig_cols), mean, na.rm = TRUE), .groups = "drop")

  cat("\nSignature scores by SPP1 group:\n")
  print(sig_comparison)

  # Transpose for easier reading
  sig_comparison_long <- sig_comparison %>%
    pivot_longer(-spp1_group, names_to = "signature", values_to = "score") %>%
    pivot_wider(names_from = spp1_group, values_from = score) %>%
    mutate(
      diff = SPP1_Pos - SPP1_Neg,
      fold_change = SPP1_Pos / SPP1_Neg
    ) %>%
    arrange(desc(diff))

  write.csv(sig_comparison_long, file.path(OUT_TABLES_SNRNASEQ, "spp1_microglia", "signature_comparison.csv"), row.names = FALSE)

  cat("\nSignature enrichment in SPP1_Pos vs SPP1_Neg:\n")
  print(sig_comparison_long)

  # Statistical test for each signature
  cat("\nStatistical comparison (Wilcoxon):\n")
  stats_list <- list()
  for (sig in sig_cols) {
    high_vals <- seu_mg@meta.data[[sig]][seu_mg$spp1_group == "SPP1_Pos"]
    low_vals <- seu_mg@meta.data[[sig]][seu_mg$spp1_group == "SPP1_Neg"]

    test <- wilcox.test(high_vals, low_vals)
    stats_list[[sig]] <- data.frame(
      signature = sig,
      spp1_high_mean = mean(high_vals, na.rm = TRUE),
      spp1_low_mean = mean(low_vals, na.rm = TRUE),
      diff = mean(high_vals, na.rm = TRUE) - mean(low_vals, na.rm = TRUE),
      pvalue = test$p.value
    )
  }

  stats_df <- bind_rows(stats_list) %>%
    mutate(fdr = p.adjust(pvalue, method = "BH")) %>%
    arrange(pvalue)

  write.csv(stats_df, file.path(OUT_TABLES_SNRNASEQ, "spp1_microglia", "signature_statistics.csv"), row.names = FALSE)

  cat("\nSignatures significantly enriched in SPP1_Pos (FDR < 0.05):\n")
  print(stats_df %>% filter(fdr < 0.05))

  # Visualization
  plot_df <- seu_mg@meta.data %>%
    filter(spp1_group %in% c("SPP1_Pos", "SPP1_Neg")) %>%
    dplyr::select(spp1_group, all_of(sig_cols)) %>%
    pivot_longer(-spp1_group, names_to = "Signature", values_to = "Score")

  # Order signatures by enrichment
  sig_order <- stats_df$signature[order(stats_df$diff, decreasing = TRUE)]
  plot_df$Signature <- factor(plot_df$Signature, levels = sig_order)

  p_sigs <- ggplot(plot_df, aes(x = Signature, y = Score, fill = spp1_group)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    scale_fill_manual(values = c("SPP1_Pos" = "#E41A1C", "SPP1_Neg" = "#377EB8")) +
    labs(title = "Microglial State Signatures: SPP1-High vs SPP1-Low",
         x = NULL, y = "Module Score", fill = "SPP1 Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = "bottom")

  ggsave(file.path(OUT_FIGURES_SNRNASEQ, "spp1_microglia", "signature_comparison_boxplot.png"),
         p_sigs, width = 14, height = 8)
}

# ============================================================
# PART 6: Condition-specific SPP1+ analysis
# ============================================================

cat("\n=== PART 6: Condition-Specific Analysis ===\n")

# Is SPP1+ phenotype more common in Implant?
spp1_cond_stats <- seu_mg@meta.data %>%
  group_by(Condition, spp1_group) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  filter(spp1_group == "SPP1_Pos")

cat("\nSPP1_Pos enrichment by condition:\n")
print(spp1_cond_stats)

# Chi-square test
contingency <- table(seu_mg$Condition, seu_mg$spp1_group)
chi_test <- chisq.test(contingency)
cat(sprintf("\nChi-square test (Condition x SPP1 group): p = %.2e\n", chi_test$p.value))

write.csv(spp1_cond_stats, file.path(OUT_TABLES_SNRNASEQ, "spp1_microglia", "spp1_high_by_condition.csv"), row.names = FALSE)

# ============================================================
# SUMMARY
# ============================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SPP1+ MICROGLIA CHARACTERIZATION SUMMARY\n")
cat(strrep("=", 60), "\n")

cat("\n[1] SPP1 Expression Distribution:\n")
cat(sprintf("    Total microglia: %d\n", ncol(seu_mg)))
cat(sprintf("    SPP1+ cells (>0): %.1f%%\n", 100 * mean(spp1_expr > 0)))
cat(sprintf("    SPP1_Pos (top 25%%): %d cells\n", sum(seu_mg$spp1_group == "SPP1_Pos")))

cat("\n[2] Differential Expression (SPP1_Pos vs SPP1_Neg):\n")
cat(sprintf("    Significant genes (FDR < 0.05): %d\n", sum(markers$p_val_adj < 0.05)))
cat(sprintf("    Upregulated: %d\n", sum(markers$p_val_adj < 0.05 & markers$avg_log2FC > 0)))
cat(sprintf("    Downregulated: %d\n", sum(markers$p_val_adj < 0.05 & markers$avg_log2FC < 0)))

if (exists("ego_up_bp") && !is.null(ego_up_bp) && nrow(ego_up_bp) > 0) {
  cat("\n[3] Top GO BP Terms (Upregulated in SPP1_Pos):\n")
  top_terms <- head(as.data.frame(ego_up_bp)$Description, 10)
  for (i in seq_along(top_terms)) {
    cat(sprintf("    %d. %s\n", i, top_terms[i]))
  }
}

if (exists("stats_df")) {
  cat("\n[4] Enriched Signatures in SPP1_Pos:\n")
  enriched <- stats_df %>% filter(fdr < 0.05, diff > 0)
  for (i in 1:min(nrow(enriched), 10)) {
    cat(sprintf("    %s: +%.3f (p = %.2e)\n",
                enriched$signature[i], enriched$diff[i], enriched$fdr[i]))
  }
}

cat("\nOutputs saved to:\n")
cat(sprintf("  Tables: %s/spp1_microglia/\n", OUT_TABLES_SNRNASEQ))
cat(sprintf("  Figures: %s/spp1_microglia/\n", OUT_FIGURES_SNRNASEQ))
