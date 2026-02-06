# 02d_concordant_timecourse.R
# Track the 50 validated concordant genes across both timecourse datasets
# Creates merged temporal trajectory of the implant signature

library(dplyr)
library(tidyr)
library(ggplot2)
source("R/config.R")

cat("=== Concordant Gene Timecourse Analysis ===\n")

# Create output directory
dir.create(file.path(OUT_FIGURES, "timecourse"), recursive = TRUE, showWarnings = FALSE)

# Load concordant genes
concordant <- read.csv(file.path(OUT_TABLES_COMPARISON, "validated_genes.csv"))
cat(sprintf("Concordant genes: %d (%d up, %d down)\n",
            nrow(concordant),
            sum(concordant$category == "Implant_Up"),
            sum(concordant$category == "Implant_Down")))

concordant_genes <- concordant$gene_upper

# ============================================================
# 1. SILICON TIMECOURSE DATA
# ============================================================

cat("\n1. Loading silicon timecourse results...\n")

# Silicon files are named de_24h_100um.csv, de_1wk_100um.csv, etc.
silicon_files <- list.files(file.path(OUT_TABLES, "timecourse"),
                            pattern = "^de_.*\\.csv",
                            full.names = TRUE)

if (length(silicon_files) == 0) {
  cat("No silicon timecourse files found. Running timecourse analysis first...\n")
  source("R/02b_timecourse_analysis.R")
  silicon_files <- list.files(file.path(OUT_TABLES, "timecourse"),
                              pattern = "^de_.*\\.csv",
                              full.names = TRUE)
}

silicon_tc <- lapply(silicon_files, function(f) {
  df <- read.csv(f)
  # Extract comparison from filename (de_24h_100um.csv -> 24h_100um)
  comp <- gsub("^de_|\\.csv$", "", basename(f))
  df$comparison <- comp
  df
}) %>% bind_rows()

cat(sprintf("Silicon: %d comparisons, %d total rows\n",
            length(unique(silicon_tc$comparison)), nrow(silicon_tc)))

# ============================================================
# 2. POLYIMIDE TIMECOURSE DATA
# ============================================================

cat("\n2. Loading polyimide timecourse results...\n")

polyimide_files <- list.files(file.path(OUT_TABLES, "timecourse"),
                              pattern = "polyimide_Week.*\\.csv",
                              full.names = TRUE)

if (length(polyimide_files) == 0) {
  cat("No polyimide timecourse files found. Running timecourse analysis first...\n")
  source("R/02c_polyimide_timecourse.R")
  polyimide_files <- list.files(file.path(OUT_TABLES, "timecourse"),
                                pattern = "polyimide_Week.*\\.csv",
                                full.names = TRUE)
}

polyimide_tc <- lapply(polyimide_files, function(f) {
  df <- read.csv(f)
  comp <- gsub("polyimide_|\\.csv", "", basename(f))
  df$comparison <- comp
  df
}) %>% bind_rows()

cat(sprintf("Polyimide: %d comparisons\n", length(unique(polyimide_tc$comparison))))

# ============================================================
# 3. EXTRACT CONCORDANT GENE LOG2FC OVER TIME
# ============================================================

cat("\n3. Extracting concordant gene temporal patterns...\n")

# Silicon - normalize gene names
silicon_tc$gene_upper <- toupper(silicon_tc$gene)
silicon_conc <- silicon_tc %>%
  filter(gene_upper %in% concordant_genes) %>%
  dplyr::select(gene_upper, log2FoldChange, padj, comparison) %>%
  rename(log2FC = log2FoldChange)

# Parse silicon comparisons
silicon_conc <- silicon_conc %>%
  mutate(
    timepoint = case_when(
      grepl("24h", comparison) ~ "24h",
      grepl("1wk", comparison) ~ "1wk",
      grepl("6wk", comparison) ~ "6wk",
      TRUE ~ "Unknown"
    ),
    distance = case_when(
      grepl("100um", comparison) ~ "100um",
      grepl("500um", comparison) ~ "500um",
      TRUE ~ "Unknown"
    ),
    platform = "Silicon"
  )

# Polyimide - normalize gene names
polyimide_tc$gene_upper <- toupper(polyimide_tc$gene)
polyimide_conc <- polyimide_tc %>%
  filter(gene_upper %in% concordant_genes) %>%
  dplyr::select(gene_upper, logFC, adj.P.Val, comparison) %>%
  rename(log2FC = logFC, padj = adj.P.Val)

# Parse polyimide comparisons
polyimide_conc <- polyimide_conc %>%
  mutate(
    timepoint = gsub("_vs_Ctrl", "", comparison),
    distance = "electrode",  # Polyimide doesn't have distance info
    platform = "Polyimide"
  )

cat(sprintf("Silicon concordant gene-timepoints: %d\n", nrow(silicon_conc)))
cat(sprintf("Polyimide concordant gene-timepoints: %d\n", nrow(polyimide_conc)))

# ============================================================
# 4. CALCULATE SIGNATURE SCORES OVER TIME
# ============================================================

cat("\n4. Calculating signature scores over time...\n")

# Get UP and DOWN gene lists
up_genes <- concordant %>% filter(category == "Implant_Up") %>% pull(gene_upper)
down_genes <- concordant %>% filter(category == "Implant_Down") %>% pull(gene_upper)

calc_signature_score <- function(df, up_genes, down_genes) {
  # Mean log2FC of UP genes minus mean log2FC of DOWN genes
  up_mean <- mean(df$log2FC[df$gene_upper %in% up_genes], na.rm = TRUE)
  down_mean <- mean(df$log2FC[df$gene_upper %in% down_genes], na.rm = TRUE)

  # Also calculate simple mean of all genes
  all_mean <- mean(df$log2FC, na.rm = TRUE)

  # Count significant
  n_sig <- sum(df$padj < 0.05, na.rm = TRUE)
  n_total <- nrow(df)

  data.frame(
    up_score = up_mean,
    down_score = down_mean,
    combined_score = up_mean - down_mean,
    mean_log2FC = all_mean,
    n_sig = n_sig,
    n_total = n_total,
    pct_sig = 100 * n_sig / n_total
  )
}

# Silicon signature scores
silicon_scores <- silicon_conc %>%
  filter(distance == "100um") %>%  # Focus on near-electrode
  group_by(timepoint, platform) %>%
  group_modify(~ calc_signature_score(.x, up_genes, down_genes)) %>%
  ungroup()

# Polyimide signature scores
polyimide_scores <- polyimide_conc %>%
  group_by(timepoint, platform) %>%
  group_modify(~ calc_signature_score(.x, up_genes, down_genes)) %>%
  ungroup()

# Combine
all_scores <- bind_rows(silicon_scores, polyimide_scores)

# Order timepoints
silicon_time_order <- c("24h", "1wk", "6wk")
polyimide_time_order <- c("Week0", "Week1", "Week2", "Week4", "Week18")

all_scores <- all_scores %>%
  mutate(
    time_numeric = case_when(
      timepoint == "24h" ~ 1/7,
      timepoint == "1wk" | timepoint == "Week1" ~ 1,
      timepoint == "Week0" ~ 0,
      timepoint == "Week2" ~ 2,
      timepoint == "Week4" ~ 4,
      timepoint == "6wk" ~ 6,
      timepoint == "Week18" ~ 18,
      TRUE ~ NA_real_
    )
  )

cat("\nSignature scores by timepoint:\n")
print(all_scores %>% dplyr::select(platform, timepoint, time_numeric, combined_score, pct_sig))

write.csv(all_scores, file.path(OUT_TABLES, "timecourse", "concordant_signature_scores.csv"), row.names = FALSE)

# ============================================================
# 5. INDIVIDUAL GENE TRAJECTORIES
# ============================================================

cat("\n5. Creating individual gene trajectories...\n")

# Combine all data
all_conc <- bind_rows(
  silicon_conc %>% filter(distance == "100um"),
  polyimide_conc
) %>%
  left_join(concordant %>% dplyr::select(gene_upper, category), by = "gene_upper") %>%
  mutate(
    time_numeric = case_when(
      timepoint == "24h" ~ 1/7,
      timepoint == "1wk" | timepoint == "Week1" ~ 1,
      timepoint == "Week0" ~ 0,
      timepoint == "Week2" ~ 2,
      timepoint == "Week4" ~ 4,
      timepoint == "6wk" ~ 6,
      timepoint == "Week18" ~ 18,
      TRUE ~ NA_real_
    )
  )

write.csv(all_conc, file.path(OUT_TABLES, "timecourse", "concordant_genes_timecourse.csv"), row.names = FALSE)

# ============================================================
# 6. FIGURES
# ============================================================

cat("\n6. Creating figures...\n")

# Figure 1: Signature score trajectory (combined)
p1 <- ggplot(all_scores, aes(x = time_numeric, y = combined_score, color = platform)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = COL_REF) +
  scale_color_manual(values = c("Silicon" = "#E41A1C", "Polyimide" = "#377EB8")) +
  scale_x_continuous(breaks = c(0, 1, 2, 4, 6, 18),
                     labels = c("0", "1", "2", "4", "6", "18")) +
  labs(title = "Concordant Signature Trajectory",
       subtitle = sprintf("Based on %d validated genes (48 up, 2 down)", nrow(concordant)),
       x = "Weeks post-implant",
       y = "Signature Score (UP - DOWN)",
       color = "Platform") +
  theme_publication() +
  theme(legend.position = c(0.85, 0.85))

save_figure(file.path(OUT_FIGURES, "timecourse", "concordant_signature_trajectory.png"), p1, width = 9, height = 6)

# Figure 2: UP genes mean trajectory
p2 <- ggplot(all_scores, aes(x = time_numeric, y = up_score, color = platform)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = COL_REF) +
  scale_color_manual(values = c("Silicon" = "#E41A1C", "Polyimide" = "#377EB8")) +
  scale_x_continuous(breaks = c(0, 1, 2, 4, 6, 18),
                     labels = c("0", "1", "2", "4", "6", "18")) +
  labs(title = "Upregulated Signature Trajectory",
       subtitle = sprintf("Mean log2FC of %d concordant UP genes", length(up_genes)),
       x = "Weeks post-implant",
       y = "Mean log2FC",
       color = "Platform") +
  theme_publication() +
  theme(legend.position = c(0.85, 0.85))

save_figure(file.path(OUT_FIGURES, "timecourse", "concordant_up_trajectory.png"), p2, width = 9, height = 6)

# Figure 3: Percent significant over time
p3 <- ggplot(all_scores, aes(x = time_numeric, y = pct_sig, color = platform)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  scale_color_manual(values = c("Silicon" = "#E41A1C", "Polyimide" = "#377EB8")) +
  scale_x_continuous(breaks = c(0, 1, 2, 4, 6, 18),
                     labels = c("0", "1", "2", "4", "6", "18")) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Concordant Signature Activation",
       subtitle = sprintf("Percent of %d genes with FDR < 0.05", nrow(concordant)),
       x = "Weeks post-implant",
       y = "% Significant",
       color = "Platform") +
  theme_publication() +
  theme(legend.position = c(0.85, 0.85))

save_figure(file.path(OUT_FIGURES, "timecourse", "concordant_pct_sig_trajectory.png"), p3, width = 9, height = 6)

# Figure 4: Key gene heatmap over time (TIME-ALIGNED)
key_genes <- c("C1QA", "C1QB", "C1QC", "C3", "GFAP", "SPP1", "TYROBP", "TREM2",
               "NPAS4", "ARC", "FOS", "EGR1")

heatmap_data <- all_conc %>%
  filter(gene_upper %in% key_genes) %>%
  mutate(gene_upper = factor(gene_upper, levels = key_genes))

# Create time-aligned labels with platform indicator
heatmap_data <- heatmap_data %>%
  mutate(
    # Create nice time labels
    time_weeks = case_when(
      time_numeric < 1 ~ sprintf("%.0fh", time_numeric * 7 * 24),
      time_numeric == 1 ~ "1wk",
      TRUE ~ sprintf("%.0fwk", time_numeric)
    ),
    # Combined label: time + platform abbreviation
    time_label = paste0(time_weeks, "\n(", substr(platform, 1, 3), ")")
  )

# Order by time_numeric first, then platform
time_order <- heatmap_data %>%
  distinct(time_label, time_numeric, platform) %>%
  arrange(time_numeric, platform) %>%
  pull(time_label)

heatmap_data$time_label <- factor(heatmap_data$time_label, levels = time_order)

p4 <- ggplot(heatmap_data, aes(x = time_label, y = gene_upper, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "log2FC",
                       limits = c(-3, 3), oob = scales::squish) +
  labs(title = "Key Concordant Genes: Time-Aligned Dynamics",
       subtitle = "Ordered by time post-implant; Sil = Silicon, Pol = Polyimide",
       x = "Time post-implant", y = NULL) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7),
        panel.grid = element_blank())

save_figure(file.path(OUT_FIGURES, "timecourse", "concordant_key_genes_heatmap.png"), p4, width = 14, height = 8)

# Figure 5: Faceted trajectory by category
category_trajectory <- all_conc %>%
  group_by(platform, timepoint, time_numeric, category) %>%
  summarize(
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    se = sd(log2FC, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

p5 <- ggplot(category_trajectory, aes(x = time_numeric, y = mean_log2FC, color = platform)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_log2FC - se, ymax = mean_log2FC + se), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = COL_REF) +
  facet_wrap(~category, scales = "free_y") +
  scale_color_manual(values = c("Silicon" = "#E41A1C", "Polyimide" = "#377EB8")) +
  scale_x_continuous(breaks = c(0, 1, 2, 4, 6, 18),
                     labels = c("0", "1", "2", "4", "6", "18")) +
  labs(title = "Concordant Signature by Category",
       subtitle = "UP genes (n=48) vs DOWN genes (n=2)",
       x = "Weeks post-implant",
       y = "Mean log2FC (± SE)",
       color = "Platform") +
  theme_publication() +
  theme(legend.position = "bottom")

save_figure(file.path(OUT_FIGURES, "timecourse", "concordant_category_trajectory.png"), p5, width = 10, height = 5)

# Figure 6: Top 10 most dynamic genes
top_dynamic <- all_conc %>%
  group_by(gene_upper) %>%
  summarize(
    range = max(log2FC, na.rm = TRUE) - min(log2FC, na.rm = TRUE),
    max_abs = max(abs(log2FC), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(range)) %>%
  head(10) %>%
  pull(gene_upper)

dynamic_data <- all_conc %>%
  filter(gene_upper %in% top_dynamic)

p6 <- ggplot(dynamic_data, aes(x = time_numeric, y = log2FC, color = gene_upper)) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = COL_REF) +
  facet_wrap(~platform, ncol = 2) +
  scale_color_brewer(palette = "Paired") +
  scale_x_continuous(breaks = c(0, 1, 2, 4, 6, 18),
                     labels = c("0", "1", "2", "4", "6", "18")) +
  labs(title = "Top 10 Most Dynamic Concordant Genes",
       subtitle = "Ranked by log2FC range across timepoints",
       x = "Weeks post-implant",
       y = "log2FC",
       color = "Gene") +
  theme_publication() +
  theme(legend.position = "right")

save_figure(file.path(OUT_FIGURES, "timecourse", "concordant_top_dynamic.png"), p6, width = 12, height = 6)

# ============================================================
# SUMMARY
# ============================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("CONCORDANT SIGNATURE TIMECOURSE SUMMARY\n")
cat(strrep("=", 60), "\n")

cat("\n[SILICON] (100μm from electrode):\n")
for (tp in silicon_time_order) {
  row <- silicon_scores %>% filter(timepoint == tp)
  if (nrow(row) > 0) {
    cat(sprintf("  %s: Score = %.2f, %d/%d (%.0f%%) significant\n",
                tp, row$combined_score, row$n_sig, row$n_total, row$pct_sig))
  }
}

cat("\n[POLYIMIDE]:\n")
for (tp in polyimide_time_order) {
  row <- polyimide_scores %>% filter(timepoint == tp)
  if (nrow(row) > 0) {
    cat(sprintf("  %s: Score = %.2f, %d/%d (%.0f%%) significant\n",
                tp, row$combined_score, row$n_sig, row$n_total, row$pct_sig))
  }
}

cat("\nSaved results to:\n")
cat(sprintf("  %s/timecourse/concordant_signature_scores.csv\n", OUT_TABLES))
cat(sprintf("  %s/timecourse/concordant_genes_timecourse.csv\n", OUT_TABLES))
cat(sprintf("  %s/timecourse/concordant_*.png (6 figures)\n", OUT_FIGURES))
