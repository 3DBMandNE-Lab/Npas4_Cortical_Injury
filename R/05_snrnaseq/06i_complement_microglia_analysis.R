# 06i_complement_microglia_analysis.R
# Deep characterization of Complement+ microglia biology
# Compare C1qa+ vs C1qa- microglia (complement-producing population)
# Parallel analysis to 06h_spp1_microglia_analysis.R

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
source("R/config.R")

cat("=== Complement+ Microglia Characterization ===\n\n")

# Create output directory
dir.create(file.path(OUT_TABLES_SNRNASEQ, "complement_microglia"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_FIGURES_SNRNASEQ, "complement_microglia"), recursive = TRUE, showWarnings = FALSE)

# Load snRNA-seq data
SNRNASEQ_SOURCE <- SNRNASEQ_PATH  # defined in config.R
cat("Loading snRNA-seq data...\n")
seu <- readRDS(SNRNASEQ_SOURCE)
cat(sprintf("Loaded: %d cells, %d genes\n", ncol(seu), nrow(seu)))

# Use RNA assay
DefaultAssay(seu) <- "RNA"

# Find microglia
celltype_col <- if ("celltype_l2" %in% colnames(seu@meta.data)) "celltype_l2" else "celltype_l1"
microglia_patterns <- c("Microglia", "microglia", "MG")
all_celltypes <- unique(seu@meta.data[[celltype_col]])
microglia_types <- all_celltypes[grepl(paste(microglia_patterns, collapse = "|"), all_celltypes, ignore.case = TRUE)]

# Subset to microglia
seu_mg <- subset(seu, cells = colnames(seu)[seu@meta.data[[celltype_col]] %in% microglia_types])
# CRITICAL: Explicitly set RNA assay after subset (subset may change default)
DefaultAssay(seu_mg) <- "RNA"
cat(sprintf("Microglia subset: %d cells (using %s assay)\n", ncol(seu_mg), DefaultAssay(seu_mg)))

# Check complement gene expression (use C1qa as primary marker)
complement_genes <- c("C1qa", "C1qb", "C1qc", "C3")
available_comp <- complement_genes[complement_genes %in% rownames(seu_mg)]
cat(sprintf("Complement genes available: %s\n", paste(available_comp, collapse = ", ")))

# Get C1qa expression (primary complement marker)
# NOTE: Use RNA assay for gene detection - SCT may alter sparse expression patterns
c1qa_expr <- GetAssayData(seu_mg, assay = "RNA", layer = "data")["C1qa", ]
cat(sprintf("C1qa expression range: %.2f - %.2f\n", min(c1qa_expr), max(c1qa_expr)))
cat(sprintf("C1qa+ cells (>0): %d (%.1f%%)\n", sum(c1qa_expr > 0), 100 * mean(c1qa_expr > 0)))

# Also get aggregate complement score
comp_genes_present <- intersect(c("C1qa", "C1qb", "C1qc", "C3"), rownames(seu_mg))
seu_mg <- AddModuleScore(seu_mg, features = list(comp_genes_present), name = "Complement_Score", seed = 42)
colnames(seu_mg@meta.data)[ncol(seu_mg@meta.data)] <- "Complement_Score"

# Define C1qa+ vs C1qa- (binary since expression is sparse)
seu_mg$c1qa_group <- ifelse(c1qa_expr > 0, "C1qa_Pos", "C1qa_Neg")
seu_mg$c1qa_expr <- c1qa_expr

# Also define complement-high based on module score (top 25%)
comp_q75 <- quantile(seu_mg$Complement_Score, 0.75)
comp_q25 <- quantile(seu_mg$Complement_Score, 0.25)
seu_mg$comp_group <- case_when(
  seu_mg$Complement_Score >= comp_q75 ~ "Comp_High",
  seu_mg$Complement_Score <= comp_q25 ~ "Comp_Low",
  TRUE ~ "Comp_Mid"
)

cat("\nC1qa group distribution:\n")
print(table(seu_mg$c1qa_group))

cat("\nComplement score group distribution:\n")
print(table(seu_mg$comp_group))

cat("\nC1qa groups by condition:\n")
print(table(seu_mg$c1qa_group, seu_mg$Condition))

cat("\nComplement score groups by condition:\n")
print(table(seu_mg$comp_group, seu_mg$Condition))

# ============================================================
# PART 1: C1qa+ enrichment by condition
# ============================================================

cat("\n=== PART 1: C1qa+ Enrichment by Condition ===\n")

c1qa_by_cond <- seu_mg@meta.data %>%
  group_by(Condition) %>%
  summarize(
    n_cells = n(),
    n_c1qa_pos = sum(c1qa_group == "C1qa_Pos"),
    pct_c1qa_pos = 100 * n_c1qa_pos / n_cells,
    mean_comp_score = mean(Complement_Score, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nC1qa expression by condition:\n")
print(c1qa_by_cond)

write.csv(c1qa_by_cond, file.path(OUT_TABLES_SNRNASEQ, "complement_microglia", "c1qa_by_condition.csv"), row.names = FALSE)

# ============================================================
# PART 2: Differential Expression - Complement High vs Low
# ============================================================

cat("\n=== PART 2: DE Analysis Complement-High vs Complement-Low ===\n")

# Use complement score groups for better power (more cells)
seu_mg_hl <- subset(seu_mg, cells = colnames(seu_mg)[seu_mg$comp_group %in% c("Comp_High", "Comp_Low")])
cat(sprintf("Cells for DE: %d Comp_High, %d Comp_Low\n",
            sum(seu_mg_hl$comp_group == "Comp_High"),
            sum(seu_mg_hl$comp_group == "Comp_Low")))

# Run Seurat FindMarkers
Idents(seu_mg_hl) <- "comp_group"
markers <- FindMarkers(seu_mg_hl, ident.1 = "Comp_High", ident.2 = "Comp_Low",
                       min.pct = 0.1, logfc.threshold = 0.25, test.use = "wilcox")

markers$gene <- rownames(markers)
markers <- markers %>% arrange(p_val_adj, desc(abs(avg_log2FC)))

cat(sprintf("\nDE results: %d genes tested\n", nrow(markers)))
cat(sprintf("Significant (FDR < 0.05): %d\n", sum(markers$p_val_adj < 0.05)))
cat(sprintf("  Upregulated in Comp_High: %d\n", sum(markers$p_val_adj < 0.05 & markers$avg_log2FC > 0)))
cat(sprintf("  Downregulated in Comp_High: %d\n", sum(markers$p_val_adj < 0.05 & markers$avg_log2FC < 0)))

write.csv(markers, file.path(OUT_TABLES_SNRNASEQ, "complement_microglia", "de_comp_high_vs_low.csv"), row.names = FALSE)

# Top genes
cat("\nTop 20 genes upregulated in Complement-High microglia:\n")
top_up <- markers %>% filter(p_val_adj < 0.05, avg_log2FC > 0) %>% head(20)
print(top_up[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")])

cat("\nTop 20 genes downregulated in Complement-High microglia:\n")
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
      cat(sprintf("\nGO BP enriched (Comp_High UP): %d terms\n", nrow(ego_up_bp)))

      ego_up_df <- as.data.frame(ego_up_bp)
      write.csv(ego_up_df, file.path(OUT_TABLES_SNRNASEQ, "complement_microglia", "go_bp_comp_high_up.csv"), row.names = FALSE)

      cat("\nTop GO BP terms (Complement-High upregulated):\n")
      print(head(ego_up_df[, c("Description", "GeneRatio", "p.adjust")], 20))

      if (nrow(ego_up_bp) >= 3) {
        p_go_up <- dotplot(ego_up_bp, showCategory = 20, title = "GO BP: Complement-High Microglia (Upregulated)")
        ggsave(file.path(OUT_FIGURES_SNRNASEQ, "complement_microglia", "go_bp_comp_high_up.png"),
               p_go_up, width = 10, height = 8)
      }
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
      cat(sprintf("\nGO BP enriched (Comp_High DOWN): %d terms\n", nrow(ego_down_bp)))

      ego_down_df <- as.data.frame(ego_down_bp)
      write.csv(ego_down_df, file.path(OUT_TABLES_SNRNASEQ, "complement_microglia", "go_bp_comp_high_down.csv"), row.names = FALSE)

      cat("\nTop GO BP terms (Complement-High downregulated):\n")
      print(head(ego_down_df[, c("Description", "GeneRatio", "p.adjust")], 15))
    }
  }
}

# ============================================================
# PART 4: Signature Analysis - Microglial States
# ============================================================

cat("\n=== PART 4: Known Microglial State Signatures ===\n")

# Same signatures as SPP1 analysis for comparison
microglial_signatures <- list(
  Homeostatic = c("P2ry12", "Tmem119", "Cx3cr1", "Siglech", "Hexb", "Cst3", "Sall1"),
  DAM_Stage1 = c("Tyrobp", "Ctsd", "Lyz2", "Ctsb", "Apoe"),
  DAM_Stage2 = c("Trem2", "Axl", "Cst7", "Lpl", "Cd9", "Csf1", "Itgax", "Clec7a", "Lilrb4"),
  Phagocytic = c("Cd68", "Mertk", "Axl", "Tyro3", "Megf10", "Gulp1"),
  Inflammatory_M1 = c("Il1b", "Tnf", "Il6", "Nos2", "Cd86", "Ccl2", "Cxcl10"),
  Anti_inflammatory_M2 = c("Arg1", "Mrc1", "Cd163", "Il10", "Tgfb1", "Chil3"),
  Proliferative = c("Mki67", "Top2a", "Pcna", "Cdk1", "Ccnb1"),
  IFN_Response = c("Ifit1", "Ifit2", "Ifit3", "Isg15", "Mx1", "Stat1", "Irf7"),
  Lipid_Metabolism = c("Apoe", "Lpl", "Fabp5", "Abca1", "Abcg1", "Plin2"),
  Synapse_Pruning = c("C1qa", "C1qb", "C1qc", "C3", "Trem2", "Tyrobp"),
  SPP1_Module = c("Spp1", "Gpnmb", "Lgals3", "Fabp5", "Igf1", "Cd63", "Lpl")
)

# Score signatures
available_genes <- rownames(seu_mg)
for (sig_name in names(microglial_signatures)) {
  genes <- microglial_signatures[[sig_name]]
  genes_title <- tools::toTitleCase(tolower(genes))
  present <- intersect(genes_title, available_genes)

  if (length(present) >= 3) {
    seu_mg <- AddModuleScore(seu_mg, features = list(present), name = sig_name, seed = 42)
    colnames(seu_mg@meta.data)[ncol(seu_mg@meta.data)] <- sig_name
  }
}

sig_cols <- intersect(names(microglial_signatures), colnames(seu_mg@meta.data))

# Compare signature scores between complement groups
if (length(sig_cols) > 0) {
  sig_comparison <- seu_mg@meta.data %>%
    filter(comp_group %in% c("Comp_High", "Comp_Low")) %>%
    group_by(comp_group) %>%
    summarize(across(all_of(sig_cols), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

  sig_comparison_long <- sig_comparison %>%
    pivot_longer(-comp_group, names_to = "signature", values_to = "score") %>%
    pivot_wider(names_from = comp_group, values_from = score) %>%
    mutate(
      diff = Comp_High - Comp_Low
    ) %>%
    arrange(desc(diff))

  write.csv(sig_comparison_long, file.path(OUT_TABLES_SNRNASEQ, "complement_microglia", "signature_comparison.csv"), row.names = FALSE)

  cat("\nSignature enrichment in Comp_High vs Comp_Low:\n")
  print(sig_comparison_long)

  # Statistical test
  stats_list <- list()
  for (sig in sig_cols) {
    high_vals <- seu_mg@meta.data[[sig]][seu_mg$comp_group == "Comp_High"]
    low_vals <- seu_mg@meta.data[[sig]][seu_mg$comp_group == "Comp_Low"]

    test <- wilcox.test(high_vals, low_vals)
    stats_list[[sig]] <- data.frame(
      signature = sig,
      comp_high_mean = mean(high_vals, na.rm = TRUE),
      comp_low_mean = mean(low_vals, na.rm = TRUE),
      diff = mean(high_vals, na.rm = TRUE) - mean(low_vals, na.rm = TRUE),
      pvalue = test$p.value
    )
  }

  stats_df <- bind_rows(stats_list) %>%
    mutate(fdr = p.adjust(pvalue, method = "BH")) %>%
    arrange(pvalue)

  write.csv(stats_df, file.path(OUT_TABLES_SNRNASEQ, "complement_microglia", "signature_statistics.csv"), row.names = FALSE)

  cat("\nSignatures significantly enriched in Comp_High (FDR < 0.05):\n")
  print(stats_df %>% filter(fdr < 0.05))

  # Visualization
  plot_df <- seu_mg@meta.data %>%
    filter(comp_group %in% c("Comp_High", "Comp_Low")) %>%
    dplyr::select(comp_group, all_of(sig_cols)) %>%
    pivot_longer(-comp_group, names_to = "Signature", values_to = "Score")

  sig_order <- stats_df$signature[order(stats_df$diff, decreasing = TRUE)]
  plot_df$Signature <- factor(plot_df$Signature, levels = sig_order)

  p_sigs <- ggplot(plot_df, aes(x = Signature, y = Score, fill = comp_group)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    scale_fill_manual(values = c("Comp_High" = "#984EA3", "Comp_Low" = "#4DAF4A")) +
    labs(title = "Microglial State Signatures: Complement-High vs Complement-Low",
         x = NULL, y = "Module Score", fill = "Complement Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = "bottom")

  ggsave(file.path(OUT_FIGURES_SNRNASEQ, "complement_microglia", "signature_comparison_boxplot.png"),
         p_sigs, width = 14, height = 8)
}

# ============================================================
# PART 5: Compare SPP1+ vs Complement+ populations
# ============================================================

cat("\n=== PART 5: SPP1+ vs Complement+ Population Overlap ===\n")

# Get SPP1 expression
spp1_expr <- GetAssayData(seu_mg, layer = "data")["Spp1", ]
seu_mg$spp1_pos <- spp1_expr > 0

# Create combined classification
seu_mg$population <- case_when(
  seu_mg$spp1_pos & seu_mg$comp_group == "Comp_High" ~ "SPP1+/Comp+",
  seu_mg$spp1_pos & seu_mg$comp_group != "Comp_High" ~ "SPP1+/Comp-",
  !seu_mg$spp1_pos & seu_mg$comp_group == "Comp_High" ~ "SPP1-/Comp+",
  TRUE ~ "SPP1-/Comp-"
)

cat("\nPopulation overlap:\n")
print(table(seu_mg$population))

cat("\nPopulation by condition:\n")
pop_by_cond <- seu_mg@meta.data %>%
  group_by(Condition, population) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(pct = 100 * n / sum(n))

print(pop_by_cond %>% arrange(Condition, population))

write.csv(pop_by_cond, file.path(OUT_TABLES_SNRNASEQ, "complement_microglia", "spp1_complement_overlap.csv"), row.names = FALSE)

# Fisher's exact test for independence
contingency <- table(seu_mg$spp1_pos, seu_mg$comp_group == "Comp_High")
fisher_test <- fisher.test(contingency)
cat(sprintf("\nFisher's exact test (SPP1+ vs Comp_High independence): p = %.4f, OR = %.2f\n",
            fisher_test$p.value, fisher_test$estimate))

# ============================================================
# SUMMARY
# ============================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("COMPLEMENT+ MICROGLIA CHARACTERIZATION SUMMARY\n")
cat(strrep("=", 60), "\n")

cat("\n[1] C1qa Expression Distribution:\n")
cat(sprintf("    Total microglia: %d\n", ncol(seu_mg)))
cat(sprintf("    C1qa+ cells (>0): %.1f%%\n", 100 * mean(c1qa_expr > 0)))

cat("\n[2] Differential Expression (Comp_High vs Comp_Low):\n")
cat(sprintf("    Significant genes (FDR < 0.05): %d\n", sum(markers$p_val_adj < 0.05)))
cat(sprintf("    Upregulated: %d\n", sum(markers$p_val_adj < 0.05 & markers$avg_log2FC > 0)))
cat(sprintf("    Downregulated: %d\n", sum(markers$p_val_adj < 0.05 & markers$avg_log2FC < 0)))

if (exists("ego_up_bp") && !is.null(ego_up_bp) && nrow(ego_up_bp) > 0) {
  cat("\n[3] Top GO BP Terms (Upregulated in Comp_High):\n")
  top_terms <- head(as.data.frame(ego_up_bp)$Description, 10)
  for (i in seq_along(top_terms)) {
    cat(sprintf("    %d. %s\n", i, top_terms[i]))
  }
}

if (exists("stats_df")) {
  cat("\n[4] Enriched Signatures in Comp_High:\n")
  enriched <- stats_df %>% filter(fdr < 0.05, diff > 0)
  for (i in 1:min(nrow(enriched), 10)) {
    cat(sprintf("    %s: +%.3f (p = %.2e)\n",
                enriched$signature[i], enriched$diff[i], enriched$fdr[i]))
  }

  cat("\n[5] Depleted Signatures in Comp_High:\n")
  depleted <- stats_df %>% filter(fdr < 0.05, diff < 0)
  for (i in 1:min(nrow(depleted), 5)) {
    cat(sprintf("    %s: %.3f (p = %.2e)\n",
                depleted$signature[i], depleted$diff[i], depleted$fdr[i]))
  }
}

cat("\n[6] SPP1+ vs Complement+ Independence:\n")
cat(sprintf("    Fisher's test p = %.4f\n", fisher_test$p.value))
cat(sprintf("    Odds ratio = %.2f\n", fisher_test$estimate))
if (fisher_test$p.value < 0.05) {
  if (fisher_test$estimate > 1) {
    cat("    Interpretation: SPP1+ and Complement+ are POSITIVELY associated\n")
  } else {
    cat("    Interpretation: SPP1+ and Complement+ are NEGATIVELY associated (distinct populations)\n")
  }
} else {
  cat("    Interpretation: SPP1+ and Complement+ are INDEPENDENT populations\n")
}

cat("\nOutputs saved to:\n")
cat(sprintf("  Tables: %s/complement_microglia/\n", OUT_TABLES_SNRNASEQ))
cat(sprintf("  Figures: %s/complement_microglia/\n", OUT_FIGURES_SNRNASEQ))
