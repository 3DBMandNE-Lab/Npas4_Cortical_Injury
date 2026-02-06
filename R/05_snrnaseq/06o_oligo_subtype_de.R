# 06o_oligo_subtype_de.R
# Pseudobulk DE for EACH oligodendrocyte subtype separately
# - Myelinating (8,791)
# - Precursor/OPC (4,972)
# - Premyelinating (659)

library(Seurat)
library(dplyr)
library(limma)
library(edgeR)
source("R/config.R")

cat("=== Oligodendrocyte Subtype-Specific DE Analysis ===\n\n")

out_dir <- file.path(OUT_TABLES_SNRNASEQ, "oligodendrocytes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
# SNRNASEQ_PATH defined in config.R
cat("Loading snRNA-seq data...\n")
seu <- readRDS(SNRNASEQ_PATH)

# CRITICAL: Use RNA assay
DefaultAssay(seu) <- "RNA"
cat(sprintf("Total cells: %d (using %s assay)\n", ncol(seu), DefaultAssay(seu)))

# Show oligodendrocyte annotations
cat("\nOligodendrocyte subtypes:\n")
print(table(seu$celltype_l2[grepl("Oligo", seu$celltype_l2)]))

# Define subtypes
subtypes <- c("Oligodendrocyte_Myelinating", "Oligodendrocyte_Precursor", "Oligodendrocyte_Premyelinating")

# Function to run pseudobulk DE
run_subtype_de <- function(seu, subtype_name, min_cells = 10) {
  cat(sprintf("\n%s\n", strrep("=", 60)))
  cat(sprintf("Analyzing: %s\n", subtype_name))
  cat(strrep("=", 60), "\n")

  cells <- subset(seu, celltype_l2 == subtype_name)
  DefaultAssay(cells) <- "RNA"

  cat(sprintf("Cells: %d\n", ncol(cells)))
  cat("By condition:\n")
  print(table(cells$Condition))

  samples <- unique(cells$orig.ident)
  counts <- GetAssayData(cells, assay = "RNA", layer = "counts")

  pb_list <- list()
  pb_meta <- data.frame()

  for (s in samples) {
    idx <- which(cells$orig.ident == s)
    if (length(idx) < min_cells) next

    pb_list[[s]] <- Matrix::rowSums(counts[, idx, drop = FALSE])
    pb_meta <- rbind(pb_meta, data.frame(
      sample = s,
      condition = cells$Condition[idx[1]],
      n_cells = length(idx)
    ))
  }

  if (length(pb_list) < 4) {
    cat("Too few samples, skipping\n")
    return(NULL)
  }

  pb_counts <- do.call(cbind, pb_list)

  impl_ctrl <- pb_meta$sample[pb_meta$condition %in% c("Implant", "Control")]
  if (length(impl_ctrl) < 4) {
    cat("Too few Implant/Control samples, skipping\n")
    return(NULL)
  }

  pb_counts <- pb_counts[, impl_ctrl]
  pb_meta <- pb_meta[pb_meta$sample %in% impl_ctrl, ]

  cat(sprintf("Samples for DE: %d (Control: %d, Implant: %d)\n",
              nrow(pb_meta),
              sum(pb_meta$condition == "Control"),
              sum(pb_meta$condition == "Implant")))

  dge <- DGEList(counts = pb_counts)
  keep <- filterByExpr(dge, group = pb_meta$condition)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)

  pb_meta$condition <- factor(pb_meta$condition, levels = c("Control", "Implant"))
  design <- model.matrix(~ condition, data = pb_meta)

  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  results <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  results$gene <- rownames(results)
  results$significant <- results$adj.P.Val < 0.05
  results$subtype <- subtype_name

  n_sig <- sum(results$significant)
  n_up <- sum(results$significant & results$logFC > 0)
  n_down <- sum(results$significant & results$logFC < 0)

  cat(sprintf("\nDE Results: %d significant (FDR < 0.05)\n", n_sig))
  cat(sprintf("  Up: %d, Down: %d\n", n_up, n_down))

  cat("\nTop 10 downregulated:\n")
  print(head(results[results$logFC < 0, c("gene", "logFC", "P.Value", "adj.P.Val")], 10))

  cat("\nTop 10 upregulated:\n")
  print(head(results[results$logFC > 0, c("gene", "logFC", "P.Value", "adj.P.Val")], 10))

  return(results)
}

# Run for each subtype
all_results <- list()

for (subtype in subtypes) {
  res <- run_subtype_de(seu, subtype)
  if (!is.null(res)) {
    all_results[[subtype]] <- res
    safe_name <- gsub("Oligodendrocyte_", "", subtype)
    write.csv(res, file.path(out_dir, paste0("de_", tolower(safe_name), "_implant_vs_ctrl.csv")),
              row.names = FALSE)
  }
}

# ============================================================
# Key oligodendrocyte genes across subtypes
# ============================================================

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("KEY OLIGODENDROCYTE GENES BY SUBTYPE\n")
cat(strrep("=", 70), "\n\n")

key_genes <- c(
  # Myelination
  "Mbp", "Plp1", "Mog", "Mag", "Cnp", "Mobp", "Olig1", "Olig2",
  # OPC markers
  "Pdgfra", "Cspg4", "Sox10",
  # MHC/Immune (from concordant signature)
  "B2m", "Cd9", "H2-D1", "H2-K1",
  # Stress/Reactive
  "Gfap", "Vim", "Serpina3n",
  # Lipid
  "Apoe", "Clu"
)

key_summary <- data.frame()
for (subtype in names(all_results)) {
  res <- all_results[[subtype]]
  key_res <- res[res$gene %in% key_genes, c("gene", "logFC", "P.Value", "adj.P.Val")]
  key_res$subtype <- gsub("Oligodendrocyte_", "", subtype)
  key_summary <- rbind(key_summary, key_res)
}

if (nrow(key_summary) > 0) {
  library(tidyr)
  wide_summary <- key_summary %>%
    select(gene, subtype, logFC) %>%
    pivot_wider(names_from = subtype, values_from = logFC)

  cat("log2FC by subtype:\n")
  print(as.data.frame(wide_summary))

  write.csv(wide_summary, file.path(out_dir, "key_genes_by_subtype.csv"), row.names = FALSE)
}

cat(sprintf("\nResults saved to: %s/\n", out_dir))
