# 06g_microglial_states.R
# Microglial state analysis from snRNA-seq
# Input: snRNA-seq reference
# Output: tables/snrnaseq/microglia_*.csv, figures/snrnaseq/microglia_*.png

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/config.R")

SNRNASEQ_SOURCE <- SNRNASEQ_PATH  # defined in config.R

cat("Loading snRNA-seq dataset...\n")
seu <- readRDS(SNRNASEQ_SOURCE)
cat(sprintf("Loaded: %d cells, %d genes\n", ncol(seu), nrow(seu)))

# Use RNA assay
DefaultAssay(seu) <- "RNA"

# Define microglial state signatures
microglial_signatures <- list(
  Homeostatic = c("P2ry12", "Tmem119", "Cx3cr1", "Siglech", "Hexb", "Cst3"),
  Activated = c("Cd68", "Cd86", "Il1b", "Tnf", "Nos2", "Ptgs2"),
  DAM = c("Trem2", "Apoe", "Lpl", "Spp1", "Tyrobp", "Cd9", "Cst7", "Itgax"),
  Inflammatory = c("Ccl2", "Ccl5", "Cxcl10", "Il6", "Cxcl1", "Ccl3")
)

# Check gene availability
available_genes <- rownames(seu)
cat("\nMicroglial signature gene availability:\n")
for (state in names(microglial_signatures)) {
  genes <- microglial_signatures[[state]]
  present <- intersect(genes, available_genes)
  cat(sprintf("  %s: %d/%d genes present\n", state, length(present), length(genes)))
}

# Subset to microglia (if celltype annotation exists)
celltype_col <- NULL
if ("celltype_l2" %in% colnames(seu@meta.data)) {
  celltype_col <- "celltype_l2"
} else if ("celltype_l1" %in% colnames(seu@meta.data)) {
  celltype_col <- "celltype_l1"
}

if (!is.null(celltype_col)) {
  cat(sprintf("\nCell types in data:\n"))
  print(table(seu@meta.data[[celltype_col]]))

  # Find microglia
  microglia_patterns <- c("Microglia", "microglia", "MG", "Immune")
  microglia_types <- unique(seu@meta.data[[celltype_col]])[
    grepl(paste(microglia_patterns, collapse = "|"), unique(seu@meta.data[[celltype_col]]), ignore.case = TRUE)
  ]

  if (length(microglia_types) > 0) {
    cat(sprintf("\nSubsetting to microglia: %s\n", paste(microglia_types, collapse = ", ")))
    seu_mg <- subset(seu, cells = colnames(seu)[seu@meta.data[[celltype_col]] %in% microglia_types])
    cat(sprintf("Microglia subset: %d cells\n", ncol(seu_mg)))
  } else {
    cat("\nNo microglia annotation found, using all cells\n")
    seu_mg <- seu
  }
} else {
  cat("\nNo cell type annotation found, using all cells\n")
  seu_mg <- seu
}

# Score microglial states
cat("\nScoring microglial states...\n")
for (state in names(microglial_signatures)) {
  genes <- microglial_signatures[[state]]
  present <- intersect(genes, available_genes)

  if (length(present) >= 2) {
    seu_mg <- AddModuleScore(seu_mg, features = list(present), name = state, seed = 42)
    colnames(seu_mg@meta.data)[ncol(seu_mg@meta.data)] <- state
    cat(sprintf("  Scored: %s (%d genes)\n", state, length(present)))
  }
}

meta <- seu_mg@meta.data
state_cols <- intersect(names(microglial_signatures), colnames(meta))

# Assign cells to dominant state
if (length(state_cols) > 0) {
  state_scores <- meta[, state_cols, drop = FALSE]

  # Assign to highest scoring state
  meta$dominant_state <- apply(state_scores, 1, function(x) {
    if (all(x < quantile(unlist(state_scores), 0.25, na.rm = TRUE))) {
      return("Low/Undefined")
    }
    names(x)[which.max(x)]
  })

  cat("\nMicroglial state distribution:\n")
  print(table(meta$dominant_state))

  # By condition
  if ("Condition" %in% colnames(meta)) {
    state_by_cond <- meta %>%
      group_by(Condition, dominant_state) %>%
      summarize(n = n(), .groups = "drop") %>%
      group_by(Condition) %>%
      mutate(pct = 100 * n / sum(n))

    cat("\nState distribution by condition:\n")
    print(state_by_cond)

    write.csv(state_by_cond, file.path(OUT_TABLES_SNRNASEQ, "microglia_state_distribution.csv"), row.names = FALSE)

    # Plot: State distribution by condition
    p1 <- ggplot(state_by_cond, aes(x = Condition, y = pct, fill = dominant_state)) +
      geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
      scale_fill_manual(values = c("Homeostatic" = "#4DAF4A", "Activated" = "#E41A1C",
                                   "DAM" = "#984EA3", "Inflammatory" = "#FF7F00",
                                   "Low/Undefined" = "grey70")) +
      labs(title = "Microglial State Distribution",
           x = NULL, y = "% of Microglia", fill = "State") +
      theme_publication()

    save_figure(file.path(OUT_FIGURES_SNRNASEQ, "microglia_state_distribution.png"), p1, width = 8, height = 6)

    # Plot: State proportions comparison
    p2 <- ggplot(state_by_cond %>% filter(dominant_state != "Low/Undefined"),
                 aes(x = dominant_state, y = pct, fill = Condition)) +
      geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
      scale_fill_manual(values = c("Control" = COL_CONTROL, "Implant" = COL_IMPLANT, "Stab" = COL_STAB)) +
      labs(title = "Microglial State Comparison",
           x = NULL, y = "% of Microglia") +
      theme_publication() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    save_figure(file.path(OUT_FIGURES_SNRNASEQ, "microglia_state_comparison.png"), p2, width = 8, height = 6)
  }

  # Summary statistics per state
  summary_list <- list()
  for (state in state_cols) {
    if ("Condition" %in% colnames(meta)) {
      state_summary <- meta %>%
        group_by(Condition) %>%
        summarize(
          mean_score = mean(.data[[state]], na.rm = TRUE),
          sd_score = sd(.data[[state]], na.rm = TRUE),
          n_cells = n(),
          .groups = "drop"
        ) %>%
        mutate(state = state)
      summary_list[[state]] <- state_summary
    }
  }

  if (length(summary_list) > 0) {
    summary_df <- bind_rows(summary_list)
    write.csv(summary_df, file.path(OUT_TABLES_SNRNASEQ, "microglia_state_scores.csv"), row.names = FALSE)

    # Statistical tests
    cat("\nStatistical comparisons (Implant vs Control):\n")
    stats_list <- list()
    for (state in state_cols) {
      ctrl_vals <- meta[[state]][meta$Condition == "Control"]
      impl_vals <- meta[[state]][meta$Condition == "Implant"]

      if (length(ctrl_vals) > 10 && length(impl_vals) > 10) {
        test <- wilcox.test(impl_vals, ctrl_vals)
        stats_list[[state]] <- data.frame(
          state = state,
          impl_mean = mean(impl_vals, na.rm = TRUE),
          ctrl_mean = mean(ctrl_vals, na.rm = TRUE),
          pvalue = test$p.value
        )
        cat(sprintf("  %s: p = %.2e (Implant %.3f vs Control %.3f)\n",
                    state, test$p.value, mean(impl_vals), mean(ctrl_vals)))
      }
    }

    if (length(stats_list) > 0) {
      stats_df <- bind_rows(stats_list)
      stats_df$fdr <- p.adjust(stats_df$pvalue, method = "BH")
      write.csv(stats_df, file.path(OUT_TABLES_SNRNASEQ, "microglia_state_statistics.csv"), row.names = FALSE)
    }

    # Violin plots of state scores
    violin_df <- meta %>%
      dplyr::select(Condition, all_of(state_cols)) %>%
      pivot_longer(-Condition, names_to = "State", values_to = "Score")

    violin_df$Condition <- factor(violin_df$Condition, levels = c("Control", "Implant", "Stab"))

    p3 <- ggplot(violin_df, aes(x = Condition, y = Score, fill = Condition)) +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
      facet_wrap(~State, scales = "free_y", ncol = 2) +
      scale_fill_manual(values = c("Control" = COL_CONTROL, "Implant" = COL_IMPLANT, "Stab" = COL_STAB)) +
      labs(title = "Microglial State Scores by Condition",
           x = NULL, y = "Module Score") +
      theme_publication() +
      theme(legend.position = "none")

    save_figure(file.path(OUT_FIGURES_SNRNASEQ, "microglia_state_violins.png"), p3, width = 8, height = 8)
  }
}

cat(sprintf("\nSaved microglial analysis to: %s/\n", OUT_TABLES_SNRNASEQ))
