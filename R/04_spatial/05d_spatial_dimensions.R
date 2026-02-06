# 05d_spatial_dimensions.R
# Quantify kill zone physical dimensions from Visium data
# Uses 10x Genomics Visium specs: 55μm spot diameter, 100μm center-to-center

library(dplyr)
library(ggplot2)
source("R/config.R")

cat("=== Spatial Domain Dimension Analysis ===\n")

# Visium specifications
SPOT_DIAMETER_UM <- 55      # μm
SPOT_SPACING_UM <- 100      # center-to-center μm
SPOT_AREA_UM2 <- pi * (SPOT_DIAMETER_UM/2)^2  # ~2376 μm²

cat(sprintf("Visium specs: %dμm spot diameter, %dμm spacing\n",
            SPOT_DIAMETER_UM, SPOT_SPACING_UM))

# Load spatial data
seurat_list <- readRDS(file.path(DATA_PROCESSED, "spatial_scored.RDS"))

results <- list()

for (sid in names(seurat_list)) {
  cat(sprintf("\n--- %s ---\n", sid))
  seu <- seurat_list[[sid]]
  meta <- seu@meta.data

  # Check for coordinates
  if (!all(c("row", "col") %in% colnames(meta))) {
    cat("  No coordinates found, skipping...\n")
    next
  }

  # Ensure coordinates are numeric
  meta$row <- as.numeric(meta$row)
  meta$col <- as.numeric(meta$col)

  # Define kill zone (high inflammatory + low neuronal)
  thresh_up <- quantile(meta$Implant_Up, 0.75, na.rm = TRUE)
  thresh_neur <- quantile(meta$Neuronal, 0.25, na.rm = TRUE)
  meta$killzone <- meta$Implant_Up > thresh_up & meta$Neuronal < thresh_neur

  n_killzone <- sum(meta$killzone, na.rm = TRUE)
  n_total <- sum(!is.na(meta$killzone))
  pct_killzone <- 100 * n_killzone / n_total

  cat(sprintf("  Kill zone spots: %d / %d (%.1f%%)\n", n_killzone, n_total, pct_killzone))

  # Calculate kill zone physical extent
  if (n_killzone > 0) {
    kz_spots <- meta[meta$killzone, ]

    # Get coordinate ranges
    row_range <- range(kz_spots$row)
    col_range <- range(kz_spots$col)

    # Convert to physical dimensions (spots * spacing)
    height_spots <- diff(row_range) + 1
    width_spots <- diff(col_range) + 1

    height_um <- height_spots * SPOT_SPACING_UM
    width_um <- width_spots * SPOT_SPACING_UM

    cat(sprintf("  Kill zone extent: %d x %d spots\n", height_spots, width_spots))
    cat(sprintf("  Kill zone extent: %d x %d μm\n", height_um, width_um))

    # Calculate approximate area
    # Method 1: Bounding box
    bbox_area_um2 <- height_um * width_um

    # Method 2: Actual spot coverage (more accurate)
    actual_area_um2 <- n_killzone * (SPOT_SPACING_UM^2)  # each spot represents 100x100 μm area

    # Method 3: Estimate radius assuming circular kill zone
    # Area = πr², so r = sqrt(Area/π)
    equiv_radius_um <- sqrt(actual_area_um2 / pi)

    cat(sprintf("  Estimated area (bbox): %.0f μm² (%.2f mm²)\n", bbox_area_um2, bbox_area_um2/1e6))
    cat(sprintf("  Estimated area (spots): %.0f μm² (%.2f mm²)\n", actual_area_um2, actual_area_um2/1e6))
    cat(sprintf("  Equivalent radius (if circular): %.0f μm\n", equiv_radius_um))

    # Calculate centroid of kill zone
    centroid_row <- mean(kz_spots$row)
    centroid_col <- mean(kz_spots$col)

    # Calculate distances from centroid
    kz_spots$dist_from_center <- sqrt((kz_spots$row - centroid_row)^2 +
                                       (kz_spots$col - centroid_col)^2) * SPOT_SPACING_UM

    # Get 90th percentile radius (most kill zone spots within this)
    radius_90pct <- quantile(kz_spots$dist_from_center, 0.90)
    radius_max <- max(kz_spots$dist_from_center)

    cat(sprintf("  90%% of kill zone within: %.0f μm radius\n", radius_90pct))
    cat(sprintf("  Max kill zone extent: %.0f μm from center\n", radius_max))

    # Compare to protected zone (low inflammatory + high neuronal)
    thresh_up_low <- quantile(meta$Implant_Up, 0.25, na.rm = TRUE)
    thresh_neur_high <- quantile(meta$Neuronal, 0.75, na.rm = TRUE)
    meta$protected <- meta$Implant_Up < thresh_up_low & meta$Neuronal > thresh_neur_high

    n_protected <- sum(meta$protected, na.rm = TRUE)
    pct_protected <- 100 * n_protected / n_total

    cat(sprintf("  Protected zone spots: %d / %d (%.1f%%)\n", n_protected, n_total, pct_protected))

    # If protected zone exists, calculate distance from kill zone centroid
    if (n_protected > 0) {
      prot_spots <- meta[meta$protected, ]
      prot_spots$dist_from_kz_center <- sqrt((prot_spots$row - centroid_row)^2 +
                                              (prot_spots$col - centroid_col)^2) * SPOT_SPACING_UM

      min_protected_dist <- min(prot_spots$dist_from_kz_center)
      median_protected_dist <- median(prot_spots$dist_from_kz_center)

      cat(sprintf("  Nearest protected zone: %.0f μm from kill zone center\n", min_protected_dist))
      cat(sprintf("  Median protected zone: %.0f μm from kill zone center\n", median_protected_dist))
    }

    # Store results
    results[[sid]] <- data.frame(
      sample = sid,
      condition = meta$condition[1],
      n_spots = n_total,
      n_killzone = n_killzone,
      killzone_pct = pct_killzone,
      killzone_height_um = height_um,
      killzone_width_um = width_um,
      killzone_area_um2 = actual_area_um2,
      killzone_equiv_radius_um = equiv_radius_um,
      killzone_90pct_radius_um = radius_90pct,
      killzone_max_radius_um = radius_max,
      n_protected = n_protected,
      protected_pct = pct_protected,
      min_protected_dist_um = if(n_protected > 0) min_protected_dist else NA,
      median_protected_dist_um = if(n_protected > 0) median_protected_dist else NA,
      stringsAsFactors = FALSE
    )
  }
}

# Combine results
results_df <- bind_rows(results)

cat("\n\n=== SUMMARY ===\n")
print(results_df[, c("sample", "condition", "killzone_equiv_radius_um",
                     "killzone_90pct_radius_um", "min_protected_dist_um")])

# Key comparison to bulk RNA-seq findings
cat("\n=== COMPARISON TO BULK RNA-SEQ KILL ZONE ===\n")
cat("Bulk RNA-seq: 100μm = catastrophic NPAS4 silencing, 500μm = protected\n")
cat(sprintf("Visium mean kill zone radius: %.0f μm (range: %.0f - %.0f)\n",
            mean(results_df$killzone_equiv_radius_um, na.rm = TRUE),
            min(results_df$killzone_equiv_radius_um, na.rm = TRUE),
            max(results_df$killzone_equiv_radius_um, na.rm = TRUE)))
cat(sprintf("Visium 90%% kill zone radius: %.0f μm (range: %.0f - %.0f)\n",
            mean(results_df$killzone_90pct_radius_um, na.rm = TRUE),
            min(results_df$killzone_90pct_radius_um, na.rm = TRUE),
            max(results_df$killzone_90pct_radius_um, na.rm = TRUE)))
cat(sprintf("Visium protected zone starts: %.0f μm from center (range: %.0f - %.0f)\n",
            mean(results_df$min_protected_dist_um, na.rm = TRUE),
            min(results_df$min_protected_dist_um, na.rm = TRUE),
            max(results_df$min_protected_dist_um, na.rm = TRUE)))

# Save results
write.csv(results_df, file.path(OUT_TABLES_SPATIAL, "killzone_dimensions.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %s/killzone_dimensions.csv\n", OUT_TABLES_SPATIAL))

# Create visualization
plot_data <- results_df %>%
  tidyr::pivot_longer(
    cols = c(killzone_equiv_radius_um, killzone_90pct_radius_um, min_protected_dist_um),
    names_to = "metric",
    values_to = "distance_um"
  ) %>%
  mutate(
    metric = factor(metric,
                    levels = c("killzone_equiv_radius_um", "killzone_90pct_radius_um", "min_protected_dist_um"),
                    labels = c("Kill Zone Radius\n(equivalent)", "Kill Zone Radius\n(90%)", "Protected Zone\nDistance"))
  )

# Add reference lines for bulk RNA-seq findings
p <- ggplot(plot_data, aes(x = condition, y = distance_um, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "#E41A1C", linewidth = 1) +
  geom_hline(yintercept = 500, linetype = "dashed", color = "#377EB8", linewidth = 1) +
  annotate("text", x = 0.5, y = 100, label = "100μm (bulk kill zone)",
           hjust = 0, vjust = -0.5, color = "#E41A1C", size = 3) +
  annotate("text", x = 0.5, y = 500, label = "500μm (bulk protected)",
           hjust = 0, vjust = -0.5, color = "#377EB8", size = 3) +
  scale_fill_manual(values = c("#E41A1C", "#FF7F00", "#377EB8")) +
  labs(title = "Spatial Kill Zone Dimensions vs Bulk RNA-seq",
       subtitle = "Visium quantification (55μm spots, 100μm spacing)",
       x = NULL, y = "Distance (μm)",
       fill = "Metric") +
  theme_publication() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

save_figure(file.path(OUT_FIGURES, "manuscript", "fig_spatial_dimensions.png"), p, width = 10, height = 7)

cat("\nSaved: output/figures/manuscript/fig_spatial_dimensions.png\n")
