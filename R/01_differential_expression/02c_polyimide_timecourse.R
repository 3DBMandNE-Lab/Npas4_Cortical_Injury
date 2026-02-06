# 02c_polyimide_timecourse.R
# Timecourse limma analysis of Polyimide microarray data
# Timepoints: Week0, Week1, Week2, Week4, Week18 (implant vs naive control)
# Input: CEL files from Polyimide_Microarray
# Output: tables/timecourse/*.csv, figures/timecourse/*.png

library(oligo)
library(limma)
library(dplyr)
library(tidyr)
library(ggplot2)
library(clariomsrathttranscriptcluster.db)
source("R/config.R")

# Create output directories
dir.create(file.path(OUT_TABLES, "timecourse"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_FIGURES, "timecourse"), recursive = TRUE, showWarnings = FALSE)

cat("=== Polyimide Timecourse Analysis ===\n\n")

# Define data paths
base_path <- "data/Internal/Polyimide_Microarray"
implant_path <- file.path(base_path, "Trepanation_Implants")
control_path <- file.path(base_path, "Naive_Controls")

# List CEL files
implant_files <- list.files(implant_path, pattern = "\\.CEL$", full.names = TRUE)
control_files <- list.files(control_path, pattern = "\\.CEL$", full.names = TRUE)

cat(sprintf("Implant files: %d\n", length(implant_files)))
cat(sprintf("Control files: %d\n", length(control_files)))

# Parse timepoints from filenames
parse_timepoint <- function(filename) {
  basename_clean <- gsub("\\.CEL$", "", basename(filename))
  if (grepl("^Week", basename_clean)) {
    week <- gsub("Week([0-9]+)_.*", "\\1", basename_clean)
    return(paste0("Week", week))
  } else if (grepl("^Control", basename_clean)) {
    return("Control")
  }
  return(NA)
}

# Create sample metadata
implant_meta <- data.frame(
  file = implant_files,
  sample = basename(implant_files),
  timepoint = sapply(implant_files, parse_timepoint),
  condition = "Implant",
  stringsAsFactors = FALSE
)

control_meta <- data.frame(
  file = control_files,
  sample = basename(control_files),
  timepoint = "Control",
  condition = "Control",
  stringsAsFactors = FALSE
)

sample_meta <- rbind(implant_meta, control_meta)
rownames(sample_meta) <- sample_meta$sample

cat("\nSample distribution:\n")
print(table(sample_meta$timepoint))

# Read CEL files
cat("\nReading CEL files...\n")
all_files <- c(implant_files, control_files)
raw_data <- read.celfiles(all_files)

# RMA normalization
cat("Running RMA normalization...\n")
eset <- rma(raw_data)
expr <- exprs(eset)

# Clean sample names
colnames(expr) <- basename(colnames(expr))

# Match sample order
sample_meta <- sample_meta[colnames(expr), ]

# Add gene annotations
cat("Adding gene annotations...\n")
probe_ids <- rownames(expr)
gene_symbols <- mapIds(clariomsrathttranscriptcluster.db,
                        keys = probe_ids,
                        column = "SYMBOL",
                        keytype = "PROBEID",
                        multiVals = "first")

# Filter to annotated probes
annotated <- !is.na(gene_symbols)
expr <- expr[annotated, ]
gene_symbols <- gene_symbols[annotated]

# Collapse to gene level (take max mean expression probe per gene)
cat("Collapsing to gene level...\n")
expr_df <- as.data.frame(expr)
expr_df$gene <- gene_symbols
expr_df$mean_expr <- rowMeans(expr_df[, colnames(expr)])

expr_gene <- expr_df %>%
  group_by(gene) %>%
  slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
  ungroup()

expr_mat <- as.matrix(expr_gene[, colnames(expr)])
rownames(expr_mat) <- expr_gene$gene

cat(sprintf("Genes after annotation: %d\n", nrow(expr_mat)))

# Set up timepoints as factors
sample_meta$timepoint <- factor(sample_meta$timepoint,
                                 levels = c("Control", "Week0", "Week1", "Week2", "Week4", "Week18"))

# limma analysis for each timepoint vs Control
cat("\nRunning limma analysis for each timepoint...\n")

design <- model.matrix(~ 0 + timepoint, data = sample_meta)
colnames(design) <- gsub("timepoint", "", colnames(design))

fit <- lmFit(expr_mat, design)

# Define contrasts: each timepoint vs Control
contrasts_list <- list(
  "Week0_vs_Ctrl" = "Week0 - Control",
  "Week1_vs_Ctrl" = "Week1 - Control",
  "Week2_vs_Ctrl" = "Week2 - Control",
  "Week4_vs_Ctrl" = "Week4 - Control",
  "Week18_vs_Ctrl" = "Week18 - Control"
)

deg_counts <- data.frame(Comparison = character(), Up = integer(), Down = integer(), Total = integer())
all_results <- list()

for (comp_name in names(contrasts_list)) {
  contrast_formula <- contrasts_list[[comp_name]]

  # Check if both groups exist
  groups_needed <- trimws(unlist(strsplit(contrast_formula, "-")))
  if (!all(groups_needed %in% colnames(design))) {
    cat(sprintf("  Skipping %s: missing groups\n", comp_name))
    next
  }

  contrast_mat <- makeContrasts(contrasts = contrast_formula, levels = design)
  fit2 <- contrasts.fit(fit, contrast_mat)
  fit2 <- eBayes(fit2)

  res <- topTable(fit2, number = Inf, adjust.method = "BH")
  res$gene <- rownames(res)

  n_up <- sum(res$adj.P.Val < FDR_THRESH & res$logFC > 0, na.rm = TRUE)
  n_down <- sum(res$adj.P.Val < FDR_THRESH & res$logFC < 0, na.rm = TRUE)

  deg_counts <- rbind(deg_counts, data.frame(
    Comparison = comp_name, Up = n_up, Down = n_down, Total = n_up + n_down
  ))

  all_results[[comp_name]] <- res
  write.csv(res, file.path(OUT_TABLES, "timecourse", paste0("polyimide_", comp_name, ".csv")), row.names = FALSE)
  cat(sprintf("  %s: %d up, %d down\n", comp_name, n_up, n_down))
}

write.csv(deg_counts, file.path(OUT_TABLES, "timecourse", "polyimide_deg_counts.csv"), row.names = FALSE)

# Key genes temporal expression
key_genes <- c("Npas4", "Arc", "Fos", "Egr1", "Nr4a1", "Nr4a3", "Bdnf",
               "C1qa", "C1qb", "C1qc", "C3", "Cfh",
               "Spp1", "Tyrobp", "Trem2", "Cd68", "Aif1", "Adgre1",
               "Gfap", "Aqp4", "Vim", "Lcn2")

key_present <- intersect(key_genes, rownames(expr_mat))
cat(sprintf("\nKey genes present: %d/%d\n", length(key_present), length(key_genes)))

if (length(key_present) > 0) {
  key_mat <- expr_mat[key_present, , drop = FALSE]

  key_expr <- as.data.frame(t(key_mat)) %>%
    tibble::rownames_to_column("Sample") %>%
    pivot_longer(-Sample, names_to = "Gene", values_to = "Expression") %>%
    left_join(sample_meta %>% dplyr::select(sample, timepoint), by = c("Sample" = "sample")) %>%
    group_by(Gene, timepoint) %>%
    dplyr::summarize(Mean = mean(Expression), SE = sd(Expression) / sqrt(n()), .groups = "drop")

  write.csv(key_expr, file.path(OUT_TABLES, "timecourse", "polyimide_key_gene_expression.csv"), row.names = FALSE)

  # Plot complement genes
  complement_genes <- c("C1qa", "C3", "Spp1")
  complement_present <- intersect(complement_genes, key_present)

  if (length(complement_present) > 0) {
    complement_data <- key_expr %>%
      filter(Gene %in% complement_present)

    p1 <- ggplot(complement_data, aes(x = timepoint, y = Mean, color = Gene, group = Gene)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
      scale_color_manual(values = c("C1qa" = "#E41A1C", "C3" = "#FF7F00", "Spp1" = "#984EA3")) +
      labs(title = "Polyimide: Complement/DAM Temporal Trajectory",
           x = NULL, y = "log2 Expression") +
      theme_publication() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    save_figure(file.path(OUT_FIGURES, "timecourse", "polyimide_complement_trajectory.png"), p1, width = 8, height = 6)
  }

  # Plot neuronal genes
  neuronal_genes <- c("Npas4", "Arc", "Fos")
  neuronal_present <- intersect(neuronal_genes, key_present)

  if (length(neuronal_present) > 0) {
    neuronal_data <- key_expr %>%
      filter(Gene %in% neuronal_present)

    p2 <- ggplot(neuronal_data, aes(x = timepoint, y = Mean, color = Gene, group = Gene)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
      scale_color_manual(values = c("Npas4" = "#377EB8", "Arc" = "#BDBDBD", "Fos" = "#4DAF4A")) +
      labs(title = "Polyimide: Neuronal IEG Temporal Trajectory",
           x = NULL, y = "log2 Expression") +
      theme_publication() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    save_figure(file.path(OUT_FIGURES, "timecourse", "polyimide_neuronal_trajectory.png"), p2, width = 8, height = 6)
  }
}

# DEG counts heatmap
if (nrow(deg_counts) > 0) {
  deg_counts$Timepoint <- gsub("_vs_Ctrl", "", deg_counts$Comparison)
  deg_counts$Timepoint <- factor(deg_counts$Timepoint, levels = c("Week0", "Week1", "Week2", "Week4", "Week18"))

  p3 <- ggplot(deg_counts, aes(x = Timepoint, y = 1, fill = Total)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%d\n(%d up)", Total, Up)), color = "white", fontface = "bold", size = 4) +
    scale_fill_gradient(low = "#FEE0D2", high = "#CB181D", name = "DEGs") +
    labs(title = "Polyimide DEG Counts Over Time",
         subtitle = "vs Naive Control, FDR < 0.05",
         x = NULL, y = NULL) +
    theme_publication() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

  save_figure(file.path(OUT_FIGURES, "timecourse", "polyimide_deg_heatmap.png"), p3, width = 10, height = 3)
}

# Compare Silicon and Polyimide timecourse overlap
silicon_1wk <- tryCatch({
  read.csv(file.path(OUT_TABLES, "timecourse", "de_1wk_100um.csv"))
}, error = function(e) NULL)

polyimide_1wk <- all_results[["Week1_vs_Ctrl"]]

if (!is.null(silicon_1wk) && !is.null(polyimide_1wk)) {
  cat("\n=== Cross-platform 1-week comparison ===\n")

  # Significant genes
  sil_sig <- silicon_1wk$gene[silicon_1wk$padj < FDR_THRESH & !is.na(silicon_1wk$padj)]
  poly_sig <- polyimide_1wk$gene[polyimide_1wk$adj.P.Val < FDR_THRESH & !is.na(polyimide_1wk$adj.P.Val)]

  overlap <- intersect(sil_sig, poly_sig)
  cat(sprintf("Silicon 1wk DEGs: %d\n", length(sil_sig)))
  cat(sprintf("Polyimide Week1 DEGs: %d\n", length(poly_sig)))
  cat(sprintf("Overlap: %d (%.1f%% of smaller set)\n",
              length(overlap), 100 * length(overlap) / min(length(sil_sig), length(poly_sig))))

  # Direction concordance
  if (length(overlap) > 0) {
    sil_dir <- setNames(sign(silicon_1wk$log2FoldChange[match(overlap, silicon_1wk$gene)]), overlap)
    poly_dir <- setNames(sign(polyimide_1wk$logFC[match(overlap, polyimide_1wk$gene)]), overlap)

    concordant <- sum(sil_dir == poly_dir, na.rm = TRUE)
    cat(sprintf("Direction concordance: %d/%d (%.1f%%)\n",
                concordant, length(overlap), 100 * concordant / length(overlap)))
  }
}

cat(sprintf("\nSaved polyimide timecourse results to: %s/timecourse/\n", OUT_TABLES))
