# 06e_snrnaseq_pseudobulk.R
# Pseudobulk differential expression by condition
# Uses the main 81,834 cell snRNA-seq dataset

library(Seurat)
library(dplyr)
library(limma)
library(edgeR)
source("R/config.R")

# Use the main snRNA-seq dataset
# SNRNASEQ_PATH defined in config.R

cat("Loading snRNA-seq data...\n")
seu <- readRDS(SNRNASEQ_PATH)
cat(sprintf("Loaded: %d cells, %d genes\n", ncol(seu), nrow(seu)))

# Use RNA assay
DefaultAssay(seu) <- "RNA"

meta <- seu@meta.data

cat("\nCondition distribution:\n")
print(table(meta$Condition))

# Get counts
counts <- GetAssayData(seu, layer = "counts")
cat(sprintf("\nData matrix: %d genes x %d cells\n", nrow(counts), ncol(counts)))

# Aggregate to pseudobulk by sample (orig.ident)
samples <- unique(meta$orig.ident)
cat(sprintf("\nSamples: %d\n", length(samples)))

# Create pseudobulk matrix
pb_list <- list()
pb_meta <- data.frame()

for (s in samples) {
  cells <- which(meta$orig.ident == s)
  if (length(cells) < 10) {
    cat(sprintf("Skipping %s (only %d cells)\n", s, length(cells)))
    next
  }

  pb_list[[s]] <- Matrix::rowSums(counts[, cells, drop = FALSE])
  pb_meta <- rbind(pb_meta, data.frame(
    sample = s,
    condition = meta$Condition[cells[1]],
    duration = meta$Duration[cells[1]],
    n_cells = length(cells)
  ))
}

pb_counts <- do.call(cbind, pb_list)
rownames(pb_meta) <- pb_meta$sample

cat("\nPseudobulk sample info:\n")
print(pb_meta)

# Filter genes with low counts
dge <- DGEList(counts = pb_counts)
keep <- filterByExpr(dge, group = pb_meta$condition)
dge <- dge[keep, , keep.lib.sizes = FALSE]
cat(sprintf("\nGenes after filtering: %d\n", nrow(dge)))

# Normalize
dge <- calcNormFactors(dge)

# Design matrix
pb_meta$condition <- factor(pb_meta$condition, levels = c("Control", "Implant", "Stab"))
design <- model.matrix(~ 0 + condition, data = pb_meta)
colnames(design) <- levels(pb_meta$condition)

# Voom transformation
v <- voom(dge, design, plot = FALSE)

# Fit model
fit <- lmFit(v, design)

# Contrasts
contrasts <- makeContrasts(
  Implant_vs_Control = Implant - Control,
  Stab_vs_Control = Stab - Control,
  Implant_vs_Stab = Implant - Stab,
  levels = design
)

fit_contrasts <- contrasts.fit(fit, contrasts)
fit_contrasts <- eBayes(fit_contrasts)

# Extract results
for (coef_name in colnames(contrasts)) {
  res <- topTable(fit_contrasts, coef = coef_name, number = Inf)
  res$gene <- rownames(res)
  res$significant <- res$adj.P.Val < FDR_THRESH & abs(res$logFC) > LFC_THRESH

  n_sig <- sum(res$significant, na.rm = TRUE)
  n_up <- sum(res$significant & res$logFC > 0, na.rm = TRUE)
  n_down <- sum(res$significant & res$logFC < 0, na.rm = TRUE)

  cat(sprintf("\n%s: %d significant (FDR < 0.05, |LFC| > 0.5) - %d up, %d down\n",
              coef_name, n_sig, n_up, n_down))

  filename <- paste0("pseudobulk_DE_", tolower(gsub("_vs_", "_vs_", coef_name)), ".csv")
  write.csv(res, file.path(OUT_TABLES_SNRNASEQ, filename), row.names = FALSE)
}

# Summary
summary_df <- data.frame(
  comparison = colnames(contrasts),
  n_genes_tested = nrow(dge),
  n_samples = nrow(pb_meta)
)

write.csv(summary_df, file.path(OUT_TABLES_SNRNASEQ, "pseudobulk_DE_summary.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %s/pseudobulk_DE_*.csv\n", OUT_TABLES_SNRNASEQ))
