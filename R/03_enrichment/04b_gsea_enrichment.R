# 04b_gsea_enrichment.R
# GSEA pathway enrichment using org.Rn.eg.db
# Input: deg/*.csv
# Output: tables/enrichment/*.csv, figures/enrichment/*.png

library(clusterProfiler)
library(org.Rn.eg.db)
library(ggplot2)
library(dplyr)
source("R/config.R")

# Create enrichment output directories
dir.create(file.path(OUT_TABLES, "enrichment"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_FIGURES, "enrichment"), recursive = TRUE, showWarnings = FALSE)

# Load DEG results
silicon <- read.csv(file.path(OUT_TABLES_DEG, "silicon_deseq2_results.csv"))
polyimide <- read.csv(file.path(OUT_TABLES_DEG, "polyimide_limma_results.csv"))

# Create ranked gene lists by log2FC (sorted descending)
# Silicon
silicon_ranked <- silicon %>%
  filter(!is.na(log2FoldChange) & !is.na(gene) & gene != "") %>%
  group_by(gene) %>%
  slice_max(abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(log2FoldChange))

silicon_genelist <- setNames(silicon_ranked$log2FoldChange, silicon_ranked$gene)
cat(sprintf("Silicon: %d unique genes for GSEA\n", length(silicon_genelist)))

# Polyimide
polyimide_ranked <- polyimide %>%
  filter(!is.na(logFC) & !is.na(gene) & gene != "") %>%
  group_by(gene) %>%
  slice_max(abs(logFC), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(logFC))

polyimide_genelist <- setNames(polyimide_ranked$logFC, polyimide_ranked$gene)
cat(sprintf("Polyimide: %d unique genes for GSEA\n", length(polyimide_genelist)))

# Run GSEA for silicon using GO terms
cat("\nRunning GSEA for silicon (GO:ALL)...\n")
gsea_silicon <- gseGO(
  geneList = silicon_genelist,
  OrgDb = org.Rn.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.25,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)
gsea_silicon_df <- as.data.frame(gsea_silicon)
cat(sprintf("Silicon: %d pathways with FDR < 0.25\n", nrow(gsea_silicon_df)))

# Run GSEA for polyimide
cat("Running GSEA for polyimide (GO:ALL)...\n")
gsea_polyimide <- gseGO(
  geneList = polyimide_genelist,
  OrgDb = org.Rn.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.25,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)
gsea_polyimide_df <- as.data.frame(gsea_polyimide)
cat(sprintf("Polyimide: %d pathways with FDR < 0.25\n", nrow(gsea_polyimide_df)))

# Save GSEA results
write.csv(gsea_silicon_df, file.path(OUT_TABLES, "enrichment", "gsea_silicon.csv"), row.names = FALSE)
write.csv(gsea_polyimide_df, file.path(OUT_TABLES, "enrichment", "gsea_polyimide.csv"), row.names = FALSE)

# Pathway-level concordance
common_pathways <- intersect(gsea_silicon_df$ID, gsea_polyimide_df$ID)
cat(sprintf("\nCommon pathways: %d\n", length(common_pathways)))

if (length(common_pathways) >= 5) {
  nes_silicon <- gsea_silicon_df$NES[match(common_pathways, gsea_silicon_df$ID)]
  nes_polyimide <- gsea_polyimide_df$NES[match(common_pathways, gsea_polyimide_df$ID)]

  nes_cor <- cor.test(nes_silicon, nes_polyimide, method = "spearman")
  cat(sprintf("NES correlation: rho = %.3f (p = %.2e)\n", nes_cor$estimate, nes_cor$p.value))

  # NES correlation plot
  nes_df <- data.frame(
    pathway = common_pathways,
    NES_silicon = nes_silicon,
    NES_polyimide = nes_polyimide,
    description = gsea_silicon_df$Description[match(common_pathways, gsea_silicon_df$ID)]
  )

  p1 <- ggplot(nes_df, aes(x = NES_silicon, y = NES_polyimide)) +
    geom_point(aes(color = NES_silicon * NES_polyimide > 0), size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = COL_REF, linetype = "dashed", se = FALSE) +
    geom_hline(yintercept = 0, color = COL_REF, linewidth = 0.5) +
    geom_vline(xintercept = 0, color = COL_REF, linewidth = 0.5) +
    scale_color_manual(values = c("FALSE" = COL_NS, "TRUE" = COL_SIG),
                       labels = c("Discordant", "Concordant")) +
    labs(title = "Pathway-Level Concordance (GO Terms)",
         subtitle = sprintf("Spearman rho = %.2f, n = %d pathways", nes_cor$estimate, length(common_pathways)),
         x = "NES (Silicon)", y = "NES (Polyimide)", color = "Direction") +
    theme_publication()

  save_figure(file.path(OUT_FIGURES, "enrichment", "gsea_nes_correlation.png"), p1, width = 8, height = 7)

  # Save concordance stats
  # NOTE: NES comparison is valid across platforms (normalized enrichment scores)
  pathway_concordance <- data.frame(
    metric = "GSEA_NES_correlation",
    rank_correlation = nes_cor$estimate,
    p_value = nes_cor$p.value,
    n_pathways = length(common_pathways)
  )
  write.csv(pathway_concordance, file.path(OUT_TABLES, "enrichment", "pathway_concordance.csv"), row.names = FALSE)
}

# Top pathways dotplot for silicon
if (nrow(gsea_silicon_df) > 0) {
  p2 <- dotplot(gsea_silicon, showCategory = 20, title = "Silicon - Enriched GO Terms", split = ".sign") +
    facet_grid(. ~ .sign) +
    theme_publication()
  save_figure(file.path(OUT_FIGURES, "enrichment", "gsea_silicon_dotplot.png"), p2, width = 12, height = 10)
}

# Top pathways dotplot for polyimide
if (nrow(gsea_polyimide_df) > 0) {
  p3 <- dotplot(gsea_polyimide, showCategory = 20, title = "Polyimide - Enriched GO Terms", split = ".sign") +
    facet_grid(. ~ .sign) +
    theme_publication()
  save_figure(file.path(OUT_FIGURES, "enrichment", "gsea_polyimide_dotplot.png"), p3, width = 12, height = 10)
}

# Combined bar plot of top pathways
top_silicon <- gsea_silicon_df %>%
  arrange(desc(abs(NES))) %>%
  head(15) %>%
  mutate(platform = "Silicon")

top_polyimide <- gsea_polyimide_df %>%
  arrange(desc(abs(NES))) %>%
  head(15) %>%
  mutate(platform = "Polyimide")

top_pathways <- bind_rows(top_silicon, top_polyimide) %>%
  mutate(direction = ifelse(NES > 0, "Up", "Down"))

if (nrow(top_pathways) > 0) {
  p4 <- ggplot(top_pathways, aes(x = reorder(Description, NES), y = NES, fill = direction)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    facet_wrap(~platform, scales = "free_y") +
    coord_flip() +
    scale_fill_manual(values = c("Up" = COL_UP, "Down" = COL_DOWN)) +
    labs(title = "Top Enriched GO Pathways",
         subtitle = "FDR < 0.25",
         x = NULL, y = "Normalized Enrichment Score") +
    theme_publication() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 8))

  save_figure(file.path(OUT_FIGURES, "enrichment", "gsea_top_pathways.png"), p4, width = 14, height = 10)
}

cat(sprintf("\nSaved enrichment results to: %s/enrichment/\n", OUT_TABLES))
