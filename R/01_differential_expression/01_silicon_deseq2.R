# 01_silicon_deseq2.R
# DESeq2 analysis: Implanted (100um, 1 week) vs Naive Control
# Output: output/tables/deg/silicon_deseq2_results.csv

library(DESeq2)
library(org.Rn.eg.db)
library(dplyr)
source("R/config.R")

# Load counts - implanted and naive control
counts_implant <- read.delim("data/Internal/Silicon_RNAseq/RS1_near_far_deseq2_raw_counts.txt", row.names = 1)
counts_naive <- read.delim("data/Internal/Silicon_RNAseq/Naivecontrols_deseq2_raw_counts.txt", row.names = 1)

# Clean column names (remove "Sample_" prefix)
colnames(counts_implant) <- gsub("Sample_", "", colnames(counts_implant))
colnames(counts_naive) <- gsub("Sample_", "", colnames(counts_naive))

# 1 Week, 100um implanted samples (from Excel metadata)
implant_samples <- c("134215", "134217", "134219", "134221", "134223", "134225")
counts_implant <- counts_implant[, implant_samples]

# All naive controls (12 samples)
naive_samples <- colnames(counts_naive)

cat(sprintf("Implanted samples (1wk, 100um): %d\n", ncol(counts_implant)))
cat(sprintf("Naive control samples: %d\n", ncol(counts_naive)))

# Merge count matrices
common_genes <- intersect(rownames(counts_implant), rownames(counts_naive))
counts_all <- cbind(counts_implant[common_genes, ], counts_naive[common_genes, ])
counts_all[is.na(counts_all)] <- 0

cat(sprintf("Total genes before filtering: %d\n", nrow(counts_all)))

# Sample metadata
sample_info <- data.frame(
  sample = colnames(counts_all),
  condition = c(rep("Implanted", length(implant_samples)),
                rep("Control", length(naive_samples))),
  row.names = colnames(counts_all)
)

cat("\nSample distribution:\n")
print(table(sample_info$condition))

# Filter: keep genes with >= 10 counts in at least 2 samples
keep <- rowSums(counts_all >= 10) >= 2
counts_all <- counts_all[keep, ]
cat(sprintf("\nGenes after filtering: %d\n", nrow(counts_all)))

# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(
  countData = counts_all,
  colData = sample_info,
  design = ~ condition
)

dds <- DESeq(dds)

# Results: Implanted vs Control
res <- results(dds, contrast = c("condition", "Implanted", "Control"), pAdjustMethod = "BH")
res_df <- as.data.frame(res)
res_df$ensembl <- rownames(res_df)

# Map Ensembl to gene symbols
anno <- AnnotationDbi::select(org.Rn.eg.db, rownames(res_df), c("SYMBOL"), "ENSEMBL")
anno <- anno[!duplicated(anno$ENSEMBL), ]
rownames(anno) <- anno$ENSEMBL
res_df$gene <- anno[res_df$ensembl, "SYMBOL"]
res_df$gene[is.na(res_df$gene)] <- res_df$ensembl[is.na(res_df$gene)]

# Significance flags
res_df$significant <- res_df$padj < FDR_THRESH & !is.na(res_df$padj)
res_df$direction <- ifelse(res_df$log2FoldChange > 0, "Up", "Down")
res_df$direction[!res_df$significant] <- "NS"

# Summary
cat(sprintf("\n=== Silicon DESeq2 Results (Implanted vs Control) ===\n"))
cat(sprintf("Significant DEGs (FDR < %.2f): %d\n", FDR_THRESH, sum(res_df$significant)))
cat(sprintf("  Up in Implant: %d\n", sum(res_df$direction == "Up")))
cat(sprintf("  Down in Implant: %d\n", sum(res_df$direction == "Down")))

# Save
write.csv(res_df, file.path(OUT_TABLES_DEG, "silicon_deseq2_results.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %s/silicon_deseq2_results.csv\n", OUT_TABLES_DEG))
