# 02b_timecourse_analysis.R
# Timecourse DESeq2 analysis of Silicon RNA-seq data (Purcell)
# Multiple timepoints: 24h, 1wk, 6wk at 100um and 500um distances
# Input: Silicon RNA-seq counts + sample metadata
# Output: tables/timecourse/*.csv, figures/timecourse/*.png

library(DESeq2)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(org.Rn.eg.db)
source("R/config.R")

# Create output directories
dir.create(file.path(OUT_TABLES, "timecourse"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_FIGURES, "timecourse"), recursive = TRUE, showWarnings = FALSE)

cat("=== Timecourse DESeq2 Analysis ===\n\n")

# Load counts
cat("Loading count data...\n")
counts_stiff <- read.delim("data/Internal/Silicon_RNAseq/RS1_near_far_deseq2_raw_counts.txt", row.names = 1)
counts_naive <- read.delim("data/Internal/Silicon_RNAseq/Naivecontrols_deseq2_raw_counts.txt", row.names = 1)

# Clean sample names
colnames(counts_stiff) <- gsub("Sample_", "", colnames(counts_stiff))
colnames(counts_naive) <- gsub("Sample_|833-EP-", "", colnames(counts_naive))

cat(sprintf("Implanted samples: %d\n", ncol(counts_stiff)))
cat(sprintf("Control samples: %d\n", ncol(counts_naive)))

# Load metadata
meta_file <- "data/Internal/Silicon_RNAseq/Purcell RNAseq Sample Identification Key.xlsx"
if (file.exists(meta_file)) {
  samples_stiff <- read_excel(meta_file, sheet = "Implanted")
  samples_naive <- read_excel(meta_file, sheet = "Control")
  samples_naive <- samples_naive[samples_naive$Timepoint == "Not applicable", ]

  # Convert IDs
  samples_stiff$`Core Sample ID` <- as.character(samples_stiff$`Core Sample ID`)
  samples_naive$`Core Sample ID` <- as.character(samples_naive$`Core Sample ID`)

  cat("\nTimepoints in implanted data:\n")
  print(table(samples_stiff$Timepoint, samples_stiff$location))

  # Subset counts to matching samples
  counts_stiff <- counts_stiff[, colnames(counts_stiff) %in% samples_stiff$`Core Sample ID`]
  counts_naive <- counts_naive[, colnames(counts_naive) %in% samples_naive$`Core Sample ID`]

  # Merge counts
  common_genes <- intersect(rownames(counts_stiff), rownames(counts_naive))
  counts_all <- cbind(counts_stiff[common_genes, ], counts_naive[common_genes, ])
  counts_all[is.na(counts_all)] <- 0

  cat(sprintf("\nMerged: %d genes, %d samples\n", nrow(counts_all), ncol(counts_all)))

  # Create metadata
  meta_stiff <- samples_stiff[samples_stiff$`Core Sample ID` %in% colnames(counts_all),
                               c("Core Sample ID", "Timepoint", "location")]
  meta_stiff$Condition <- "Implanted"

  meta_naive <- data.frame(
    `Core Sample ID` = samples_naive$`Core Sample ID`[samples_naive$`Core Sample ID` %in% colnames(counts_all)],
    Timepoint = "Control",
    location = "Control",
    Condition = "Control",
    check.names = FALSE
  )

  samples_all <- rbind(as.data.frame(meta_stiff), as.data.frame(meta_naive))
  rownames(samples_all) <- samples_all$`Core Sample ID`
  samples_all <- samples_all[colnames(counts_all), ]

  samples_all$Timepoint <- factor(samples_all$Timepoint,
                                   levels = c("Control", "24 Hours", "1 Week", "6 Weeks"))
  samples_all$Distance <- factor(samples_all$location,
                                  levels = c("Control", "100um", "500um"))
  samples_all$Group <- factor(paste(samples_all$Timepoint, samples_all$Distance, sep = "_"))

  cat("\nSample groups:\n")
  print(table(samples_all$Group))

  # Filter and run DESeq2
  keep <- rowSums(counts_all >= 10) >= 2
  counts_filtered <- counts_all[keep, ]
  cat(sprintf("\nGenes after filtering: %d\n", nrow(counts_filtered)))

  # DESeq2
  cat("\nRunning DESeq2...\n")
  dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(counts_filtered),
    colData = samples_all,
    design = ~ Group
  )
  dds <- DESeq(dds)
  norm_counts <- counts(dds, normalized = TRUE)

  # Gene symbols
  gene_symbols <- mapIds(org.Rn.eg.db, keys = rownames(norm_counts),
                         column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

  # Pairwise comparisons vs Control
  comparisons <- list(
    "24h_100um" = c("Group", "24 Hours_100um", "Control_Control"),
    "24h_500um" = c("Group", "24 Hours_500um", "Control_Control"),
    "1wk_100um" = c("Group", "1 Week_100um", "Control_Control"),
    "1wk_500um" = c("Group", "1 Week_500um", "Control_Control"),
    "6wk_100um" = c("Group", "6 Weeks_100um", "Control_Control"),
    "6wk_500um" = c("Group", "6 Weeks_500um", "Control_Control")
  )

  deg_counts <- data.frame(Comparison = character(), Up = integer(), Down = integer(), Total = integer())

  for (comp_name in names(comparisons)) {
    if (all(comparisons[[comp_name]][2:3] %in% levels(samples_all$Group))) {
      res <- results(dds, contrast = comparisons[[comp_name]], alpha = FDR_THRESH)
      res_df <- as.data.frame(res)
      res_df$gene <- gene_symbols[rownames(res_df)]
      res_df$ensembl <- rownames(res_df)

      n_up <- sum(res_df$padj < FDR_THRESH & res_df$log2FoldChange > 0, na.rm = TRUE)
      n_down <- sum(res_df$padj < FDR_THRESH & res_df$log2FoldChange < 0, na.rm = TRUE)

      deg_counts <- rbind(deg_counts, data.frame(
        Comparison = comp_name, Up = n_up, Down = n_down, Total = n_up + n_down
      ))

      write.csv(res_df, file.path(OUT_TABLES, "timecourse", paste0("de_", comp_name, ".csv")), row.names = FALSE)
      cat(sprintf("  %s: %d up, %d down\n", comp_name, n_up, n_down))
    }
  }

  write.csv(deg_counts, file.path(OUT_TABLES, "timecourse", "deg_counts_by_condition.csv"), row.names = FALSE)

  # Key gene expression for figures
  key_genes <- c("Npas4", "Arc", "Fos", "Egr1", "Nr4a1", "Nr4a3", "Bdnf",
                 "C1qa", "C1qb", "C1qc", "C3", "Cfh",
                 "Spp1", "Tyrobp", "Trem2", "Cd68", "Aif1", "Adgre1",
                 "Gfap", "Aqp4", "Vim", "Lcn2")

  key_ensembl <- names(gene_symbols)[gene_symbols %in% key_genes]
  key_ensembl <- key_ensembl[!is.na(key_ensembl)]

  if (length(key_ensembl) > 0) {
    key_mat <- norm_counts[key_ensembl, , drop = FALSE]
    rownames(key_mat) <- gene_symbols[rownames(key_mat)]

    key_expr <- as.data.frame(key_mat) %>%
      tibble::rownames_to_column("Gene") %>%
      pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
      left_join(samples_all %>% tibble::rownames_to_column("Sample") %>%
                  dplyr::select(Sample, Timepoint, Distance), by = "Sample") %>%
      group_by(Gene, Timepoint, Distance) %>%
      dplyr::summarize(Mean = mean(Expression), SE = sd(Expression) / sqrt(n()), .groups = "drop")

    write.csv(key_expr, file.path(OUT_TABLES, "timecourse", "key_gene_expression.csv"), row.names = FALSE)

    # Figures
    # 1. Temporal trajectory for complement genes
    complement_genes <- c("C1qa", "C3", "Spp1")
    complement_data <- key_expr %>%
      filter(Gene %in% complement_genes, Distance == "100um")

    if (nrow(complement_data) > 0) {
      p1 <- ggplot(complement_data, aes(x = Timepoint, y = Mean, color = Gene, group = Gene)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
        scale_color_manual(values = c("C1qa" = "#E41A1C", "C3" = "#FF7F00", "Spp1" = "#984EA3")) +
        labs(title = "Complement/DAM Temporal Trajectory",
             subtitle = "100μm from electrode",
             x = NULL, y = "Normalized Expression") +
        theme_publication()

      save_figure(file.path(OUT_FIGURES, "timecourse", "complement_trajectory.png"), p1, width = 8, height = 6)
    }

    # 2. Neuronal gene trajectory
    neuronal_genes <- c("Npas4", "Arc", "Fos")
    neuronal_data <- key_expr %>%
      filter(Gene %in% neuronal_genes, Distance == "100um")

    if (nrow(neuronal_data) > 0) {
      p2 <- ggplot(neuronal_data, aes(x = Timepoint, y = Mean, color = Gene, group = Gene)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
        scale_color_manual(values = c("Npas4" = "#377EB8", "Arc" = "#BDBDBD", "Fos" = "#4DAF4A")) +
        labs(title = "Neuronal IEG Temporal Trajectory",
             subtitle = "100μm from electrode",
             x = NULL, y = "Normalized Expression") +
        theme_publication()

      save_figure(file.path(OUT_FIGURES, "timecourse", "neuronal_trajectory.png"), p2, width = 8, height = 6)
    }

    # 3. DEG count heatmap
    if (nrow(deg_counts) > 0) {
      deg_matrix <- deg_counts %>%
        separate(Comparison, into = c("Time", "Dist"), sep = "_", remove = FALSE) %>%
        mutate(Time = factor(Time, levels = c("24h", "1wk", "6wk")),
               Dist = factor(Dist, levels = c("100um", "500um")))

      p3 <- ggplot(deg_matrix, aes(x = Time, y = Dist, fill = Total)) +
        geom_tile(color = "white", linewidth = 1) +
        geom_text(aes(label = sprintf("%d\n(%d↑)", Total, Up)), color = "white", fontface = "bold", size = 4) +
        scale_fill_gradient(low = "#FEE0D2", high = "#CB181D", name = "DEGs") +
        labs(title = "DEG Counts by Condition",
             subtitle = "FDR < 0.05",
             x = "Timepoint", y = "Distance") +
        theme_publication()

      save_figure(file.path(OUT_FIGURES, "timecourse", "deg_heatmap.png"), p3, width = 7, height = 5)
    }

    # 4. Kill zone trajectory (C1qa vs Npas4)
    if (all(c("C1qa", "Npas4") %in% key_expr$Gene)) {
      # Get control baselines
      ctrl_c1qa <- key_expr$Mean[key_expr$Gene == "C1qa" & key_expr$Distance == "Control"]
      ctrl_npas4 <- key_expr$Mean[key_expr$Gene == "Npas4" & key_expr$Distance == "Control"]

      if (length(ctrl_c1qa) > 0 && length(ctrl_npas4) > 0) {
        traj_data <- key_expr %>%
          filter(Gene %in% c("C1qa", "Npas4")) %>%
          dplyr::select(Gene, Timepoint, Distance, Mean) %>%
          pivot_wider(names_from = Gene, values_from = Mean) %>%
          mutate(
            C1qa_FC = C1qa / ctrl_c1qa[1],
            Npas4_FC = pmax(Npas4 / ctrl_npas4[1], 0.01)
          ) %>%
          filter(Distance %in% c("100um", "500um"))

        p4 <- ggplot(traj_data, aes(x = C1qa_FC, y = Npas4_FC, color = Distance)) +
          geom_path(aes(group = Distance), linewidth = 1.2,
                    arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
          geom_point(aes(size = as.numeric(Timepoint)), alpha = 0.8) +
          geom_hline(yintercept = 1, linetype = "dashed", color = COL_REF) +
          geom_vline(xintercept = 1, linetype = "dashed", color = COL_REF) +
          scale_x_log10() +
          scale_y_log10() +
          scale_color_manual(values = c("100um" = COL_UP, "500um" = COL_STAB)) +
          scale_size_continuous(range = c(2, 6), guide = "none") +
          annotate("rect", xmin = 1, xmax = Inf, ymin = 0, ymax = 1, fill = COL_UP, alpha = 0.1) +
          annotate("text", x = 10, y = 0.05, label = "Kill Zone", color = COL_UP, fontface = "bold") +
          labs(title = "Kill Zone Trajectory",
               subtitle = "Complement up + Neuronal down",
               x = "C1qa Fold Change", y = "Npas4 Fold Change") +
          theme_publication()

        save_figure(file.path(OUT_FIGURES, "timecourse", "killzone_trajectory.png"), p4, width = 8, height = 7)
      }
    }
  }
} else {
  cat("ERROR: Metadata file not found:", meta_file, "\n")
}

cat(sprintf("\nSaved timecourse results to: %s/timecourse/\n", OUT_TABLES))
