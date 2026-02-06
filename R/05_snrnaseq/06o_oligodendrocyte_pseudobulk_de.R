# 06o_oligodendrocyte_pseudobulk_de.R
# Pseudobulk differential expression for OLIGODENDROCYTES ONLY
# 14,422 cells - significant but undercharacterized population

library(Seurat)
library(dplyr)
library(limma)
library(edgeR)
library(clusterProfiler)
library(org.Rn.eg.db)
source("R/config.R")

cat("=== Oligodendrocyte Pseudobulk DE Analysis ===\n\n")

# Create output directory
out_dir <- file.path(OUT_TABLES_SNRNASEQ, "oligodendrocytes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
# SNRNASEQ_PATH defined in config.R
cat("Loading snRNA-seq data...\n")
seu <- readRDS(SNRNASEQ_PATH)

# CRITICAL: Use RNA assay for all analyses
DefaultAssay(seu) <- "RNA"
cat(sprintf("Total cells: %d (using %s assay)\n", ncol(seu), DefaultAssay(seu)))

# ============================================================
# Subset to oligodendrocytes only
# ============================================================

cat("\nSubsetting to oligodendrocytes...\n")

# Check what cell type annotations are available
cat("Cell type L1 distribution:\n")
print(table(seu$celltype_l1))

# Subset to all oligodendrocyte lineage cells
oligo <- subset(seu, celltype_l1 %in% c("Oligodendrocyte", "OPC",
                                         "Oligodendrocyte_Myelinating",
                                         "Oligodendrocyte_Precursor",
                                         "Oligodendrocyte_Premyelinating"))

# If celltype_l1 doesn't work, try celltype_l2
if (ncol(oligo) < 1000) {
  cat("Trying celltype_l2...\n")
  print(table(seu$celltype_l2))
  oligo <- subset(seu, grepl("Oligo|OPC|Myelin", celltype_l2, ignore.case = TRUE))
}

# Ensure RNA assay after subset
DefaultAssay(oligo) <- "RNA"
cat(sprintf("\nOligodendrocytes: %d cells (using %s assay)\n", ncol(oligo), DefaultAssay(oligo)))

cat("\nOligodendrocytes by condition:\n")
print(table(oligo$Condition))

cat("\nOligodendrocyte subtypes:\n")
if ("celltype_l2" %in% colnames(oligo@meta.data)) {
  print(table(oligo$celltype_l2))
}

# ============================================================
# Pseudobulk aggregation by sample
# ============================================================

cat("\n=== Creating oligodendrocyte pseudobulk ===\n")

samples <- unique(oligo$orig.ident)
cat(sprintf("Number of samples: %d\n", length(samples)))

# Get counts from RNA assay explicitly
counts <- GetAssayData(oligo, assay = "RNA", layer = "counts")

pseudobulk_list <- list()
sample_meta <- data.frame()

for (s in samples) {
  cells <- colnames(oligo)[oligo$orig.ident == s]
  if (length(cells) < 20) {  # Need enough oligos per sample
    cat(sprintf("  Skipping %s (only %d oligodendrocytes)\n", s, length(cells)))
    next
  }

  pseudobulk_list[[s]] <- Matrix::rowSums(counts[, cells, drop = FALSE])

  condition <- unique(oligo$Condition[oligo$orig.ident == s])
  sample_meta <- rbind(sample_meta, data.frame(
    sample = s,
    condition = condition[1],
    n_cells = length(cells)
  ))
}

pseudobulk_matrix <- do.call(cbind, pseudobulk_list)
colnames(pseudobulk_matrix) <- names(pseudobulk_list)

cat(sprintf("\nPseudobulk matrix: %d genes x %d samples\n",
            nrow(pseudobulk_matrix), ncol(pseudobulk_matrix)))
cat("Samples by condition:\n")
print(table(sample_meta$condition))

# ============================================================
# limma-voom: Implant vs Control (oligodendrocytes only)
# ============================================================

cat("\n=== Running limma-voom: Implant vs Control (oligodendrocytes only) ===\n")

# Filter to Implant and Control
impl_ctrl_samples <- sample_meta$sample[sample_meta$condition %in% c("Implant", "Control")]
impl_ctrl_matrix <- pseudobulk_matrix[, impl_ctrl_samples]
impl_ctrl_meta <- sample_meta[sample_meta$condition %in% c("Implant", "Control"), ]

cat(sprintf("Implant vs Control: %d samples\n", nrow(impl_ctrl_meta)))
print(table(impl_ctrl_meta$condition))

# Create DGEList and filter
dge <- DGEList(counts = impl_ctrl_matrix)
keep <- filterByExpr(dge, group = impl_ctrl_meta$condition)
dge <- dge[keep, , keep.lib.sizes = FALSE]
cat(sprintf("Genes after filtering: %d\n", nrow(dge)))

# Normalize and voom
dge <- calcNormFactors(dge)
impl_ctrl_meta$condition <- factor(impl_ctrl_meta$condition, levels = c("Control", "Implant"))
design <- model.matrix(~ condition, data = impl_ctrl_meta)
colnames(design) <- c("Intercept", "Implant_vs_Control")

v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Get results
results <- topTable(fit, coef = "Implant_vs_Control", number = Inf, sort.by = "P")
results$gene <- rownames(results)
results$significant <- results$adj.P.Val < 0.05

n_sig <- sum(results$significant, na.rm = TRUE)
n_up <- sum(results$significant & results$logFC > 0, na.rm = TRUE)
n_down <- sum(results$significant & results$logFC < 0, na.rm = TRUE)

cat(sprintf("\nOligodendrocyte DE results:\n"))
cat(sprintf("  Significant (FDR < 0.05): %d genes\n", n_sig))
cat(sprintf("  Upregulated: %d\n", n_up))
cat(sprintf("  Downregulated: %d\n", n_down))

# Save full results
write.csv(results, file.path(out_dir, "de_implant_vs_control.csv"), row.names = FALSE)

# ============================================================
# Show top DE genes
# ============================================================

cat("\n=== Top Downregulated Genes in Oligodendrocytes (Implant) ===\n")
top_down <- results %>%
  filter(logFC < 0) %>%
  arrange(P.Value) %>%
  head(20)

print(top_down[, c("gene", "logFC", "P.Value", "adj.P.Val")])

cat("\n=== Top Upregulated Genes in Oligodendrocytes (Implant) ===\n")
top_up <- results %>%
  filter(logFC > 0) %>%
  arrange(P.Value) %>%
  head(20)

print(top_up[, c("gene", "logFC", "P.Value", "adj.P.Val")])

# ============================================================
# Check key oligodendrocyte genes
# ============================================================

cat("\n=== Key Oligodendrocyte Genes ===\n")
key_genes <- c(
  # Myelination
  "Mbp", "Plp1", "Mog", "Mag", "Cnp", "Mobp", "Olig1", "Olig2",
  # OPC markers
  "Pdgfra", "Cspg4", "Sox10",
  # MHC/Immune (concordant signature)
  "B2m", "Cd9", "H2-D1", "H2-K1",
  # Stress response
  "Hspa1a", "Hspa1b", "Atf3",
  # Lipid metabolism
  "Apoe", "Clu", "Lpl"
)

key_results <- results[results$gene %in% key_genes, ]
key_results <- key_results[order(key_results$logFC), ]

cat("\nKey genes sorted by log2FC:\n")
print(key_results[, c("gene", "logFC", "P.Value", "adj.P.Val", "significant")])

write.csv(key_results, file.path(out_dir, "key_oligo_genes.csv"), row.names = FALSE)

# ============================================================
# GO Enrichment (if enough significant genes)
# ============================================================

if (n_down > 20) {
  cat("\n=== GO Enrichment for Downregulated Genes ===\n")

  down_genes <- results$gene[results$significant & results$logFC < 0]

  gene_map <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Rn.eg.db)

  if (nrow(gene_map) > 10) {
    ego_down <- enrichGO(gene = gene_map$ENTREZID,
                         OrgDb = org.Rn.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2)

    if (!is.null(ego_down) && nrow(ego_down) > 0) {
      cat("\nTop GO terms (downregulated in oligodendrocytes):\n")
      print(head(ego_down@result[, c("Description", "GeneRatio", "p.adjust", "Count")], 15))
      write.csv(ego_down@result, file.path(out_dir, "go_bp_down.csv"), row.names = FALSE)
    }
  }
}

if (n_up > 20) {
  cat("\n=== GO Enrichment for Upregulated Genes ===\n")

  up_genes <- results$gene[results$significant & results$logFC > 0]

  gene_map <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Rn.eg.db)

  if (nrow(gene_map) > 10) {
    ego_up <- enrichGO(gene = gene_map$ENTREZID,
                       OrgDb = org.Rn.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

    if (!is.null(ego_up) && nrow(ego_up) > 0) {
      cat("\nTop GO terms (upregulated in oligodendrocytes):\n")
      print(head(ego_up@result[, c("Description", "GeneRatio", "p.adjust", "Count")], 15))
      write.csv(ego_up@result, file.path(out_dir, "go_bp_up.csv"), row.names = FALSE)
    }
  }
}

# ============================================================
# Summary
# ============================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("OLIGODENDROCYTE PSEUDOBULK DE SUMMARY\n")
cat(strrep("=", 60), "\n\n")

cat(sprintf("Total oligodendrocytes analyzed: %d\n", ncol(oligo)))
cat(sprintf("Samples used: %d\n", nrow(impl_ctrl_meta)))
cat(sprintf("Genes tested: %d\n", nrow(results)))
cat(sprintf("Significant (FDR < 0.05): %d\n", n_sig))
cat(sprintf("  - Downregulated: %d\n", n_down))
cat(sprintf("  - Upregulated: %d\n", n_up))

cat(sprintf("\nOutput saved to: %s/\n", out_dir))
