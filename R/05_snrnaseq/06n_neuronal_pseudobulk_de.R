# 06n_neuronal_pseudobulk_de.R
# Pseudobulk differential expression for NEURONS ONLY
# This is the proper analysis we're missing - neuron-specific DE
# Mirrors the approach used for SPP1+ and Complement+ microglia

library(Seurat)
library(dplyr)
library(limma)
library(edgeR)
library(clusterProfiler)
library(org.Rn.eg.db)
source("R/config.R")

cat("=== Neuronal Pseudobulk DE Analysis ===\n\n")

# Create output directory
out_dir <- file.path(OUT_TABLES_SNRNASEQ, "neurons")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
# SNRNASEQ_PATH defined in config.R
cat("Loading snRNA-seq data...\n")
seu <- readRDS(SNRNASEQ_PATH)

# CRITICAL: Use RNA assay for all analyses
DefaultAssay(seu) <- "RNA"
cat(sprintf("Total cells: %d (using %s assay)\n", ncol(seu), DefaultAssay(seu)))

# ============================================================
# Subset to neurons only
# ============================================================

cat("\nSubsetting to neurons...\n")

# Check what cell type annotations are available
cat("Cell type L1 distribution:\n")
print(table(seu$celltype_l1))

# Subset to neurons (excitatory + inhibitory)
neurons <- subset(seu, celltype_l1 %in% c("Neuron_Excitatory", "Neuron_Inhibitory",
                                           "Excitatory", "Inhibitory",
                                           "Neuron"))

# If celltype_l1 doesn't work, try celltype_l2
if (ncol(neurons) < 1000) {
  cat("Trying celltype_l2...\n")
  print(table(seu$celltype_l2))
  neurons <- subset(seu, grepl("Neuron|Excit|Inhib|PV|CCK|SST|VIP", celltype_l2, ignore.case = TRUE))
}

# Ensure RNA assay after subset
DefaultAssay(neurons) <- "RNA"
cat(sprintf("\nNeurons: %d cells (using %s assay)\n", ncol(neurons), DefaultAssay(neurons)))

cat("\nNeurons by condition:\n")
print(table(neurons$Condition))

cat("\nNeurons by subtype:\n")
if ("celltype_l2" %in% colnames(neurons@meta.data)) {
  print(table(neurons$celltype_l2))
}

# ============================================================
# Pseudobulk aggregation by sample
# ============================================================

cat("\n=== Creating neuronal pseudobulk ===\n")

samples <- unique(neurons$orig.ident)
cat(sprintf("Number of samples: %d\n", length(samples)))

# Get counts from RNA assay explicitly
counts <- GetAssayData(neurons, assay = "RNA", layer = "counts")

pseudobulk_list <- list()
sample_meta <- data.frame()

for (s in samples) {
  cells <- colnames(neurons)[neurons$orig.ident == s]
  if (length(cells) < 50) {  # Need enough neurons per sample
    cat(sprintf("  Skipping %s (only %d neurons)\n", s, length(cells)))
    next
  }

  pseudobulk_list[[s]] <- Matrix::rowSums(counts[, cells, drop = FALSE])

  condition <- unique(neurons$Condition[neurons$orig.ident == s])
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
# limma-voom: Implant vs Control (neurons only)
# ============================================================

cat("\n=== Running limma-voom: Implant vs Control (neurons only) ===\n")

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

cat(sprintf("\nNeuronal DE results:\n"))
cat(sprintf("  Significant (FDR < 0.05): %d genes\n", n_sig))
cat(sprintf("  Upregulated: %d\n", n_up))
cat(sprintf("  Downregulated: %d\n", n_down))

# Save full results
write.csv(results, file.path(out_dir, "de_implant_vs_control.csv"), row.names = FALSE)

# ============================================================
# Show top DE genes
# ============================================================

cat("\n=== Top Downregulated Genes in Neurons (Implant) ===\n")
top_down <- results %>%
  filter(logFC < 0) %>%
  arrange(P.Value) %>%
  head(20)

print(top_down[, c("gene", "logFC", "P.Value", "adj.P.Val")])

cat("\n=== Top Upregulated Genes in Neurons (Implant) ===\n")
top_up <- results %>%
  filter(logFC > 0) %>%
  arrange(P.Value) %>%
  head(20)

print(top_up[, c("gene", "logFC", "P.Value", "adj.P.Val")])

# ============================================================
# Check key neuronal genes
# ============================================================

cat("\n=== Key Neuronal Genes ===\n")
key_genes <- c("Npas4", "Arc", "Fos", "Egr1", "Nr4a1",  # IEGs
               "Rbfox3", "Gria1", "Gria2", "Shank2", "Homer1",  # Synaptic
               "Camk2a", "Camk2b", "Syt1", "Snap25",  # Activity
               "Nefl", "Nefm", "Nefh", "Map2",  # Structural
               "Bdnf", "Ntrk2",  # Survival
               "Slc17a7", "Grin1", "Grin2a",  # Glutamatergic
               "Gad1", "Gad2", "Pvalb", "Sst", "Cck")  # GABAergic

key_results <- results[results$gene %in% key_genes, ]
key_results <- key_results[order(key_results$logFC), ]

cat("\nKey genes sorted by log2FC:\n")
print(key_results[, c("gene", "logFC", "P.Value", "adj.P.Val", "significant")])

write.csv(key_results, file.path(out_dir, "key_neuronal_genes.csv"), row.names = FALSE)

# ============================================================
# GO Enrichment (if enough significant genes)
# ============================================================

if (n_down > 20) {
  cat("\n=== GO Enrichment for Downregulated Genes ===\n")

  down_genes <- results$gene[results$significant & results$logFC < 0]

  # Convert to Entrez IDs
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
      cat("\nTop GO terms (downregulated in neurons):\n")
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
      cat("\nTop GO terms (upregulated in neurons):\n")
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
cat("NEURONAL PSEUDOBULK DE SUMMARY\n")
cat(strrep("=", 60), "\n\n")

cat(sprintf("Total neurons analyzed: %d\n", ncol(neurons)))
cat(sprintf("Samples used: %d\n", nrow(impl_ctrl_meta)))
cat(sprintf("Genes tested: %d\n", nrow(results)))
cat(sprintf("Significant (FDR < 0.05): %d\n", n_sig))
cat(sprintf("  - Downregulated: %d\n", n_down))
cat(sprintf("  - Upregulated: %d\n", n_up))

cat(sprintf("\nOutput saved to: %s/\n", out_dir))
