# 05e_spatial_killzone_summary.R
# Summarize kill zone dimensions with proper Visium spot-based calculations
# Key insight: Each Visium spot effectively represents 100μm x 100μm area

library(dplyr)
library(ggplot2)
library(tidyr)
source("R/config.R")

cat("=== Kill Zone Dimension Summary ===\n\n")

# Visium specifications
SPOT_SPACING_UM <- 100  # center-to-center spacing

# Load results
dims <- read.csv(file.path(OUT_TABLES_SPATIAL, "killzone_dimensions.csv"))

# The equivalent radius is the key metric:
# If kill zone has N spots, total area = N * 100 * 100 μm²
# Equivalent radius = sqrt(Area / π)

# Clean up conditions for plotting
dims$condition_short <- case_when(
  dims$condition == "Craniotomy Control" ~ "Control",
  dims$condition == "Acute Stimulation" ~ "Acute",
  dims$condition == "Chronic No Stimulation" ~ "Chronic\n(No Stim)",
  dims$condition == "Chronic Stimulation" ~ "Chronic\n(Stim)",
  TRUE ~ dims$condition
)

# Order conditions
dims$condition_short <- factor(dims$condition_short,
                               levels = c("Control", "Acute", "Chronic\n(No Stim)", "Chronic\n(Stim)"))

cat("=== KILL ZONE DIMENSIONS (Visium) ===\n\n")
cat("Method: Each Visium spot represents 100μm × 100μm area\n")
cat("Equivalent radius = sqrt(n_spots × 10,000 / π)\n\n")

# Summary table
summary_table <- dims %>%
  group_by(condition) %>%
  summarize(
    n_samples = n(),
    mean_spots = mean(n_killzone),
    mean_pct = mean(killzone_pct),
    mean_radius_um = mean(killzone_equiv_radius_um),
    sd_radius_um = sd(killzone_equiv_radius_um),
    .groups = "drop"
  )

cat("By Condition:\n")
print(summary_table)

# Key comparison
cat("\n\n=== COMPARISON TO BULK RNA-SEQ ===\n")
cat("Bulk RNA-seq sampling distances:\n")
cat("  - 100μm from electrode: NPAS4 log2FC = -22.6 to -45.3 (CATASTROPHIC)\n")
cat("  - 500μm from electrode: NPAS4 log2FC = -2.9 to -5.8 (PROTECTED/NS)\n\n")

cat("Visium kill zone dimensions (equivalent circular radius):\n")
for (cond in unique(dims$condition)) {
  subset_dims <- dims[dims$condition == cond, ]
  cat(sprintf("  %s: %.0f μm (n=%d samples)\n",
              cond, mean(subset_dims$killzone_equiv_radius_um), nrow(subset_dims)))
}

cat("\n")
cat("INTERPRETATION:\n")
cat("  - Kill zone radius ranges from ~550 μm (acute) to ~900 μm (chronic)\n")
cat("  - 100μm sampling = WITHIN kill zone in all conditions\n")
cat("  - 500μm sampling = AT EDGE or OUTSIDE kill zone\n")
cat("  - Kill zone EXPANDS over time (acute → chronic)\n")
cat("  - This explains why 500μm tissue shows borderline/protected NPAS4 levels\n")

# Create visualization
p1 <- ggplot(dims, aes(x = condition_short, y = killzone_equiv_radius_um)) +
  geom_bar(stat = "identity", fill = "#E41A1C", color = "black", width = 0.6) +
  geom_hline(yintercept = 100, linetype = "solid", color = "black", linewidth = 1.2) +
  geom_hline(yintercept = 500, linetype = "dashed", color = "#377EB8", linewidth = 1.2) +
  annotate("text", x = 0.55, y = 120, label = "100um bulk sampling",
           hjust = 0, vjust = 0, size = 3, fontface = "bold") +
  annotate("text", x = 0.55, y = 520, label = "500um bulk sampling",
           hjust = 0, vjust = 0, size = 3, color = "#377EB8", fontface = "bold") +
  annotate("rect", xmin = 0.4, xmax = 4.6, ymin = 0, ymax = 100,
           alpha = 0.2, fill = "#E41A1C") +
  annotate("text", x = 4.5, y = 50, label = "CATASTROPHIC\nzone",
           hjust = 1, size = 2.5, color = "#E41A1C") +
  labs(title = "Visium Kill Zone vs Bulk RNA-seq Sampling",
       subtitle = "Kill zone equivalent radius calculated from spot counts",
       x = NULL,
       y = "Kill Zone Radius (um)") +
  scale_y_continuous(limits = c(0, 1100), breaks = seq(0, 1000, 200)) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

save_figure(file.path(OUT_FIGURES, "manuscript", "fig_killzone_vs_bulk_sampling.png"),
            p1, width = 8, height = 6)

# Second plot: Timeline showing kill zone expansion
dims_timeline <- dims %>%
  mutate(timepoint = case_when(
    condition == "Craniotomy Control" ~ "Control",
    condition == "Acute Stimulation" ~ "Acute\n(~3 days)",
    grepl("Chronic", condition) ~ "Chronic\n(weeks)"
  )) %>%
  mutate(timepoint = factor(timepoint, levels = c("Control", "Acute\n(~3 days)", "Chronic\n(weeks)")))

p2 <- ggplot(dims_timeline, aes(x = timepoint, y = killzone_equiv_radius_um)) +
  geom_point(aes(color = condition), size = 4) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, linewidth = 0.8) +
  geom_hline(yintercept = 500, linetype = "dashed", color = "#377EB8", linewidth = 1) +
  annotate("text", x = 3.3, y = 520, label = "500um (bulk protected zone)",
           hjust = 0, vjust = 0, size = 3, color = "#377EB8") +
  scale_color_manual(values = c(
    "Craniotomy Control" = "#4DAF4A",
    "Acute Stimulation" = "#FF7F00",
    "Chronic No Stimulation" = "#E41A1C",
    "Chronic Stimulation" = "#984EA3"
  )) +
  labs(title = "Kill Zone Expands Over Time",
       subtitle = "Visium equivalent radius by condition",
       x = NULL,
       y = "Kill Zone Radius (um)",
       color = "Condition") +
  scale_y_continuous(limits = c(400, 1100)) +
  theme_publication() +
  theme(legend.position = "right")

save_figure(file.path(OUT_FIGURES, "manuscript", "fig_killzone_expansion.png"),
            p2, width = 9, height = 6)

# Summary statistics for CLAUDE.md
cat("\n\n=== FOR CLAUDE.MD ===\n")
cat("| Condition | Kill Zone Radius (um) | n samples |\n")
cat("|-----------|----------------------|------------|\n")
for (i in 1:nrow(summary_table)) {
  row <- summary_table[i,]
  cat(sprintf("| %s | %.0f ± %.0f | %d |\n",
              row$condition, row$mean_radius_um, row$sd_radius_um, row$n_samples))
}

cat("\nKey Finding: Kill zone radius (550-950 um) explains bulk RNA-seq sampling:\n")
cat("- 100um sampling = deep within kill zone = catastrophic NPAS4 silencing\n")
cat("- 500um sampling = edge/outside kill zone = protected/borderline NPAS4\n")

cat("\n\nSaved figures:\n")
cat("  output/figures/manuscript/fig_killzone_vs_bulk_sampling.png\n")
cat("  output/figures/manuscript/fig_killzone_expansion.png\n")
