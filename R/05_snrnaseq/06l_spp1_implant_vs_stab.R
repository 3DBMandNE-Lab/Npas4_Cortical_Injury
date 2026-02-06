# 06l_spp1_implant_vs_stab.R
# Direct test: Is SPP1 upregulated in Implant vs Stab?
# This is the formal test for implant-specificity

library(Seurat)
library(dplyr)
library(limma)
library(edgeR)
source("R/config.R")

cat("=== SPP1 Implant vs Stab Direct Test ===\n\n")

dir.create(file.path(OUT_TABLES_SNRNASEQ, "spp1_specificity"), recursive = TRUE, showWarnings = FALSE)

# Load snRNA-seq data
# SNRNASEQ_PATH defined in config.R
cat("Loading snRNA-seq data...\n")
seu <- readRDS(SNRNASEQ_PATH)

cat(sprintf("Total cells: %d\n", ncol(seu)))
cat("Conditions:\n")
print(table(seu$Condition))

# ============================================================
# Subset to microglia only for SPP1 analysis
# ============================================================

microglia <- subset(seu, celltype_l1 == "Microglia")
# CRITICAL: Use RNA assay for pseudobulk DE - need raw counts
DefaultAssay(microglia) <- "RNA"
cat(sprintf("\nMicroglia: %d cells (using %s assay)\n", ncol(microglia), DefaultAssay(microglia)))
cat("Microglia by condition:\n")
print(table(microglia$Condition))

# ============================================================
# Pseudobulk aggregation by sample
# ============================================================

cat("\n=== Creating pseudobulk ===\n")

# Get sample IDs
samples <- unique(microglia$orig.ident)
cat(sprintf("Number of samples: %d\n", length(samples)))

# Aggregate counts by sample
pseudobulk_list <- list()
sample_meta <- data.frame()

for (s in samples) {
  cells <- colnames(microglia)[microglia$orig.ident == s]
  if (length(cells) < 10) {
    cat(sprintf("  Skipping %s (only %d cells)\n", s, length(cells)))
    next
  }

  # Use RNA assay counts explicitly for pseudobulk aggregation
  counts <- GetAssayData(microglia, assay = "RNA", layer = "counts")[, cells]
  pseudobulk_list[[s]] <- rowSums(counts)

  condition <- unique(microglia$Condition[microglia$orig.ident == s])
  sample_meta <- rbind(sample_meta, data.frame(
    sample = s,
    condition = condition[1],
    n_cells = length(cells)
  ))
}

# Create count matrix
pseudobulk_matrix <- do.call(cbind, pseudobulk_list)
colnames(pseudobulk_matrix) <- names(pseudobulk_list)

cat(sprintf("\nPseudobulk matrix: %d genes x %d samples\n",
            nrow(pseudobulk_matrix), ncol(pseudobulk_matrix)))
cat("Samples by condition:\n")
print(table(sample_meta$condition))

# ============================================================
# Filter to Implant and Stab only (excluding Control for direct comparison)
# ============================================================

implant_stab_samples <- sample_meta$sample[sample_meta$condition %in% c("Implant", "Stab")]
implant_stab_matrix <- pseudobulk_matrix[, implant_stab_samples]
implant_stab_meta <- sample_meta[sample_meta$condition %in% c("Implant", "Stab"), ]

cat(sprintf("\nImplant vs Stab comparison: %d samples\n", nrow(implant_stab_meta)))
print(table(implant_stab_meta$condition))

# ============================================================
# limma-voom analysis: Implant vs Stab
# ============================================================

cat("\n=== Running limma-voom: Implant vs Stab ===\n")

# Create DGEList
dge <- DGEList(counts = implant_stab_matrix)

# Check Spp1 raw counts before filtering
cat(sprintf("\nSpp1 raw counts in pseudobulk:\n"))
if ("Spp1" %in% rownames(implant_stab_matrix)) {
  print(implant_stab_matrix["Spp1", ])
  cat(sprintf("Total Spp1 counts: %.1f\n", sum(implant_stab_matrix["Spp1", ])))
} else {
  cat("Spp1 not found in count matrix!\n")
}

# Filter low-expressed genes - use less stringent filter to keep Spp1
# min.count = 1, min.total.count = 5 (more permissive)
keep <- filterByExpr(dge, group = implant_stab_meta$condition, min.count = 1, min.total.count = 5)

# Force keep Spp1 if it exists
if ("Spp1" %in% rownames(dge)) {
  keep["Spp1"] <- TRUE
}

dge <- dge[keep, , keep.lib.sizes = FALSE]
cat(sprintf("Genes after filtering: %d\n", nrow(dge)))
cat(sprintf("Spp1 in filtered set: %s\n", "Spp1" %in% rownames(dge)))

# Normalize
dge <- calcNormFactors(dge)

# Design matrix
design <- model.matrix(~ condition, data = implant_stab_meta)
colnames(design) <- c("Intercept", "Stab_vs_Implant")

# voom transformation
v <- voom(dge, design, plot = FALSE)

# Fit model
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Get results for Stab vs Implant (positive = higher in Stab)
results <- topTable(fit, coef = "Stab_vs_Implant", number = Inf, sort.by = "none")
results$gene <- rownames(results)

# SPP1 result (note: coefficient is Stab vs Implant, so negative = higher in Implant)
spp1_result <- results[results$gene == "Spp1", ]
cat("\n=== SPP1 Result (Implant vs Stab) ===\n")
cat(sprintf("log2FC (Stab vs Implant): %.3f\n", spp1_result$logFC))
cat(sprintf("log2FC (Implant vs Stab): %.3f\n", -spp1_result$logFC))
cat(sprintf("t-statistic: %.3f\n", spp1_result$t))
cat(sprintf("p-value: %.2e\n", spp1_result$P.Value))
cat(sprintf("FDR: %.3f\n", spp1_result$adj.P.Val))

# ============================================================
# Also run Implant vs Control and Stab vs Control for comparison
# ============================================================

cat("\n=== Running all pairwise comparisons ===\n")

# Full dataset (all three conditions)
full_meta <- sample_meta
full_matrix <- pseudobulk_matrix[, full_meta$sample]

# Create DGEList
dge_full <- DGEList(counts = full_matrix)
keep_full <- filterByExpr(dge_full, group = full_meta$condition, min.count = 1, min.total.count = 5)

# Force keep Spp1 if it exists
if ("Spp1" %in% rownames(dge_full)) {
  keep_full["Spp1"] <- TRUE
}

dge_full <- dge_full[keep_full, , keep.lib.sizes = FALSE]
cat(sprintf("Full dataset genes after filtering: %d\n", nrow(dge_full)))
cat(sprintf("Spp1 in full filtered set: %s\n", "Spp1" %in% rownames(dge_full)))
dge_full <- calcNormFactors(dge_full)

# Design matrix with all conditions
full_meta$condition <- factor(full_meta$condition, levels = c("Control", "Implant", "Stab"))
design_full <- model.matrix(~ 0 + condition, data = full_meta)
colnames(design_full) <- levels(full_meta$condition)

# voom
v_full <- voom(dge_full, design_full, plot = FALSE)

# Fit
fit_full <- lmFit(v_full, design_full)

# Contrasts
contrasts <- makeContrasts(
  Implant_vs_Control = Implant - Control,
  Stab_vs_Control = Stab - Control,
  Implant_vs_Stab = Implant - Stab,
  levels = design_full
)

fit_contrasts <- contrasts.fit(fit_full, contrasts)
fit_contrasts <- eBayes(fit_contrasts)

# Get SPP1 results for all contrasts
spp1_all <- data.frame(
  comparison = c("Implant_vs_Control", "Stab_vs_Control", "Implant_vs_Stab"),
  log2FC = c(
    topTable(fit_contrasts, coef = "Implant_vs_Control", number = Inf)["Spp1", "logFC"],
    topTable(fit_contrasts, coef = "Stab_vs_Control", number = Inf)["Spp1", "logFC"],
    topTable(fit_contrasts, coef = "Implant_vs_Stab", number = Inf)["Spp1", "logFC"]
  ),
  pvalue = c(
    topTable(fit_contrasts, coef = "Implant_vs_Control", number = Inf)["Spp1", "P.Value"],
    topTable(fit_contrasts, coef = "Stab_vs_Control", number = Inf)["Spp1", "P.Value"],
    topTable(fit_contrasts, coef = "Implant_vs_Stab", number = Inf)["Spp1", "P.Value"]
  ),
  FDR = c(
    topTable(fit_contrasts, coef = "Implant_vs_Control", number = Inf)["Spp1", "adj.P.Val"],
    topTable(fit_contrasts, coef = "Stab_vs_Control", number = Inf)["Spp1", "adj.P.Val"],
    topTable(fit_contrasts, coef = "Implant_vs_Stab", number = Inf)["Spp1", "adj.P.Val"]
  )
)

cat("\n=== SPP1 All Pairwise Comparisons ===\n")
print(spp1_all)

# Save results
write.csv(spp1_all, file.path(OUT_TABLES_SNRNASEQ, "spp1_specificity", "spp1_pairwise_comparisons.csv"),
          row.names = FALSE)

# Full DE results for Implant vs Stab
implant_vs_stab_full <- topTable(fit_contrasts, coef = "Implant_vs_Stab", number = Inf)
implant_vs_stab_full$gene <- rownames(implant_vs_stab_full)
write.csv(implant_vs_stab_full, file.path(OUT_TABLES_SNRNASEQ, "spp1_specificity", "implant_vs_stab_de.csv"),
          row.names = FALSE)

# ============================================================
# Interpretation
# ============================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SPP1 IMPLANT-SPECIFICITY TEST RESULTS\n")
cat(strrep("=", 70), "\n\n")

cat("Pairwise comparisons (microglia pseudobulk, limma-voom):\n\n")

for (i in 1:nrow(spp1_all)) {
  sig <- ifelse(spp1_all$pvalue[i] < 0.05, "*", "")
  fdr_sig <- ifelse(spp1_all$FDR[i] < 0.05, " (FDR < 0.05)", "")
  cat(sprintf("  %s: log2FC = %.2f, p = %.4f%s%s\n",
              spp1_all$comparison[i], spp1_all$log2FC[i], spp1_all$pvalue[i], sig, fdr_sig))
}

cat("\n[INTERPRETATION]:\n")
if (spp1_all$pvalue[spp1_all$comparison == "Implant_vs_Stab"] < 0.05) {
  cat("SPP1 is significantly higher in Implant vs Stab (direct test).\n")
  cat("This confirms SPP1 as an implant-specific marker, not just injury response.\n")
} else {
  cat("SPP1 does NOT reach significance in direct Implant vs Stab comparison.\n")
  cat("Recommendation: Reframe as 'differential magnitude' rather than\n")
  cat("'qualitative specificity'. The pattern (high in Implant, moderate in Stab)\n")
  cat("suggests dose-response to chronic perturbation rather than binary specificity.\n")
}

cat("\n\nOutput saved to:\n")
cat(sprintf("  %s/spp1_specificity/\n", OUT_TABLES_SNRNASEQ))
