# 04c_cell_deconvolution.R
# Cell type deconvolution of bulk RNA-seq using snRNA-seq reference
# Input: DEG results, snRNA-seq reference
# Output: tables/deconvolution/*.csv, figures/deconvolution/*.png

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/config.R")

# Create output directories
dir.create(file.path(OUT_TABLES, "deconvolution"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_FIGURES, "deconvolution"), recursive = TRUE, showWarnings = FALSE)

# Load snRNA-seq reference
SNRNASEQ_SOURCE <- SNRNASEQ_PATH  # defined in config.R
cat("Loading snRNA-seq reference...\n")
seu <- readRDS(SNRNASEQ_SOURCE)
cat(sprintf("Reference: %d cells, %d genes\n", ncol(seu), nrow(seu)))

# Use RNA assay
DefaultAssay(seu) <- "RNA"

# Get cell type column
celltype_col <- if ("celltype_l2" %in% colnames(seu@meta.data)) "celltype_l2" else "seurat_clusters"
cat(sprintf("Using cell type column: %s\n", celltype_col))
print(table(seu@meta.data[[celltype_col]]))

# Get marker genes for each cell type
cat("\nFinding cell type markers...\n")
Idents(seu) <- celltype_col
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,
                          max.cells.per.ident = 1000, verbose = FALSE)
cat(sprintf("Found %d marker genes\n", nrow(markers)))

# Top 50 markers per cell type
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50) %>%
  ungroup()

write.csv(top_markers, file.path(OUT_TABLES, "deconvolution", "celltype_markers.csv"), row.names = FALSE)

# Create signature matrix
cat("Creating signature matrix...\n")
celltypes <- unique(seu@meta.data[[celltype_col]])
expr_data <- GetAssayData(seu, layer = "data")

# Load concordant genes
cat("\nAnalyzing concordant gene cell type attribution...\n")

# Try validated_genes.csv first, fall back to shared_degs.csv
concordant_file <- file.path(OUT_TABLES, "comparison", "validated_genes.csv")
if (!file.exists(concordant_file)) {
  concordant_file <- file.path(OUT_TABLES, "comparison", "shared_degs.csv")
}
concordant <- read.csv(concordant_file)
cat(sprintf("Concordant genes: %d\n", nrow(concordant)))

# Get gene column
if ("gene_impl" %in% colnames(concordant)) {
  concordant_genes <- as.character(concordant$gene_impl)
} else if ("gene" %in% colnames(concordant)) {
  concordant_genes <- as.character(concordant$gene)
} else {
  concordant_genes <- as.character(concordant[[1]])
}

# Match to snRNA-seq reference
concordant_in_ref <- intersect(concordant_genes, rownames(expr_data))
cat(sprintf("Concordant genes in reference: %d\n", length(concordant_in_ref)))

if (length(concordant_in_ref) > 0) {
  # Calculate mean expression per cell type for concordant genes
  attribution_results <- data.frame()

  for (gene in concordant_in_ref) {
    gene_expr <- as.numeric(expr_data[gene, ])

    ct_expr <- sapply(celltypes, function(ct) {
      cells_idx <- which(seu@meta.data[[celltype_col]] == ct)
      if (length(cells_idx) > 0) {
        mean(gene_expr[cells_idx], na.rm = TRUE)
      } else {
        0
      }
    })

    # Handle edge cases
    total_expr <- sum(ct_expr, na.rm = TRUE)
    if (is.na(total_expr) || total_expr == 0) next

    ct_pct <- ct_expr / total_expr * 100
    ct_pct[is.na(ct_pct)] <- 0

    primary_ct <- names(ct_pct)[which.max(ct_pct)]
    if (length(primary_ct) == 0) next

    row_data <- data.frame(
      gene = gene,
      primary_celltype = primary_ct,
      primary_pct = max(ct_pct, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    # Add cell type percentages
    for (ct in names(ct_pct)) {
      row_data[[gsub(" ", "_", ct)]] <- ct_pct[ct]
    }

    attribution_results <- bind_rows(attribution_results, row_data)
  }

  cat(sprintf("Successfully attributed: %d genes\n", nrow(attribution_results)))

  if (nrow(attribution_results) > 0) {
    write.csv(attribution_results,
              file.path(OUT_TABLES, "deconvolution", "concordant_celltype_attribution.csv"),
              row.names = FALSE)

    cat("\nCell type attribution summary:\n")
    print(table(attribution_results$primary_celltype))

    # Plot: Cell type attribution pie chart
    ct_counts <- as.data.frame(table(attribution_results$primary_celltype))
    colnames(ct_counts) <- c("CellType", "Count")
    ct_counts$Pct <- ct_counts$Count / sum(ct_counts$Count) * 100

    # Clean up cell type names for plotting
    ct_counts$Label <- gsub("_", "\n", ct_counts$CellType)

    p1 <- ggplot(ct_counts, aes(x = "", y = Count, fill = CellType)) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar("y") +
      geom_text(aes(label = sprintf("%.1f%%", Pct)),
                position = position_stack(vjust = 0.5), size = 3, color = "white") +
      labs(title = "Concordant Gene Cell Type Attribution",
           subtitle = sprintf("n = %d genes", nrow(attribution_results))) +
      theme_void() +
      theme(legend.title = element_blank())

    save_figure(file.path(OUT_FIGURES, "deconvolution", "concordant_attribution_pie.png"), p1, width = 8, height = 8)

    # Classify genes by neuronal vs glial
    neuronal_types <- c("Neuron_Excitatory", "Neuron_Inhibitory")
    glial_types <- c("Astrocyte", "Microglia", "Oligodendrocyte_Myelinating",
                     "Oligodendrocyte_Precursor", "Oligodendrocyte_Premyelinating")

    attribution_results$category <- case_when(
      attribution_results$primary_celltype %in% neuronal_types ~ "Neuronal",
      attribution_results$primary_celltype %in% glial_types ~ "Glial",
      TRUE ~ "Other"
    )

    cat("\nNeuronal vs Glial attribution:\n")
    print(table(attribution_results$category))

    # Bar plot of neuronal vs glial
    cat_summary <- as.data.frame(table(attribution_results$category))
    colnames(cat_summary) <- c("Category", "Count")

    p2 <- ggplot(cat_summary, aes(x = Category, y = Count, fill = Category)) +
      geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.4) +
      geom_text(aes(label = Count), vjust = -0.3, size = 4) +
      scale_fill_manual(values = c("Neuronal" = COL_DOWN, "Glial" = COL_UP, "Other" = COL_NS)) +
      labs(title = "Concordant Genes: Neuronal vs Glial Origin",
           subtitle = sprintf("n = %d genes", nrow(attribution_results)),
           x = NULL, y = "Number of Genes") +
      theme_publication() +
      theme(legend.position = "none")

    save_figure(file.path(OUT_FIGURES, "deconvolution", "concordant_neuronal_glial.png"), p2, width = 6, height = 5)
  }
}

cat(sprintf("\nSaved deconvolution results to: %s/deconvolution/\n", OUT_TABLES))
