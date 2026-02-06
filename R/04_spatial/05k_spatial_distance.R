# 05k_spatial_distance.R
# Distance-dependent spatial analysis
# Input: data/processed/spatial/*.RDS
# Output: tables/spatial/distance_*.csv, figures/spatial/distance_*.png

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/config.R")

# Load spatial data
seurat_list <- readRDS(file.path(DATA_PROCESSED, "spatial", "seurat_list.RDS"))
cat(sprintf("Loaded %d spatial samples\n", length(seurat_list)))

# Define signatures to analyze
signatures_to_analyze <- c("Implant_Up", "DAM", "Complement", "Neuronal", "Astrocyte")

# Function to estimate electrode center from inflammation scores
find_electrode_center <- function(seu) {
  meta <- seu@meta.data
  coords <- GetTissueCoordinates(seu)

  # Use high inflammation as proxy for electrode location
  if ("Implant_Up" %in% colnames(meta)) {
    score_col <- "Implant_Up"
  } else if ("DAM" %in% colnames(meta)) {
    score_col <- "DAM"
  } else {
    return(NULL)
  }

  # Find centroid of top 5% inflammation spots
  thresh <- quantile(meta[[score_col]], 0.95, na.rm = TRUE)
  high_inflam <- meta[[score_col]] >= thresh

  if (sum(high_inflam) < 5) {
    return(NULL)
  }

  center_x <- mean(coords$x[high_inflam], na.rm = TRUE)
  center_y <- mean(coords$y[high_inflam], na.rm = TRUE)

  return(c(x = center_x, y = center_y))
}

# Function to calculate distance from electrode
calc_distance_from_electrode <- function(seu, center) {
  coords <- GetTissueCoordinates(seu)
  distance <- sqrt((coords$x - center["x"])^2 + (coords$y - center["y"])^2)
  return(distance)
}

# Process each sample
results_list <- list()
all_distance_data <- list()

for (sample_name in names(seurat_list)) {
  seu <- seurat_list[[sample_name]]
  meta <- seu@meta.data

  cat(sprintf("\nProcessing %s (%d spots)...\n", sample_name, ncol(seu)))

  # Find electrode center
  center <- find_electrode_center(seu)

  if (is.null(center)) {
    cat("  Could not identify electrode center, skipping\n")
    next
  }

  cat(sprintf("  Electrode center: (%.1f, %.1f)\n", center["x"], center["y"]))

  # Calculate distance from electrode
  distance <- calc_distance_from_electrode(seu, center)

  # Convert to approximate microns (Visium spots are ~100um apart)
  # Spot spacing is typically 100um center-to-center
  pixel_to_um <- 100 / median(diff(sort(unique(GetTissueCoordinates(seu)$x))))
  distance_um <- distance * pixel_to_um

  # Bin by distance
  meta$distance_um <- distance_um
  meta$distance_bin <- cut(distance_um,
                           breaks = c(0, 100, 200, 500, 1000, Inf),
                           labels = c("0-100", "100-200", "200-500", "500-1000", ">1000"),
                           include.lowest = TRUE)

  # Get available signatures
  avail_sigs <- intersect(signatures_to_analyze, colnames(meta))

  if (length(avail_sigs) == 0) {
    cat("  No signatures found, skipping\n")
    next
  }

  # Calculate mean score per distance bin
  for (sig in avail_sigs) {
    bin_summary <- meta %>%
      group_by(distance_bin) %>%
      summarize(
        mean_score = mean(.data[[sig]], na.rm = TRUE),
        sd_score = sd(.data[[sig]], na.rm = TRUE),
        n_spots = n(),
        .groups = "drop"
      ) %>%
      mutate(
        sample = sample_name,
        signature = sig,
        condition = meta$condition[1]
      )

    results_list[[paste(sample_name, sig, sep = "_")]] <- bin_summary
  }

  # Store individual spot data for detailed analysis
  spot_data <- meta %>%
    dplyr::select(distance_um, distance_bin, all_of(avail_sigs)) %>%
    mutate(sample = sample_name, condition = meta$condition[1])

  all_distance_data[[sample_name]] <- spot_data
}

# Combine results
if (length(results_list) > 0) {
  distance_summary <- bind_rows(results_list)
  write.csv(distance_summary, file.path(OUT_TABLES_SPATIAL, "distance_signature_summary.csv"), row.names = FALSE)

  cat(sprintf("\n\nDistance analysis summary:\n"))
  print(distance_summary %>% filter(signature == "Implant_Up") %>% arrange(sample, distance_bin))

  # Plot: Signature scores by distance
  p1 <- ggplot(distance_summary, aes(x = distance_bin, y = mean_score, color = signature, group = signature)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_score - sd_score/sqrt(n_spots),
                      ymax = mean_score + sd_score/sqrt(n_spots)),
                  width = 0.2) +
    facet_wrap(~sample, scales = "free_y") +
    scale_color_manual(values = COL_SPATIAL) +
    labs(title = "Signature Scores by Distance from Electrode",
         x = "Distance (μm)", y = "Mean Score") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_figure(file.path(OUT_FIGURES_SPATIAL, "distance_signatures_by_sample.png"), p1, width = 14, height = 10)

  # Average across samples by condition
  condition_summary <- distance_summary %>%
    group_by(condition, distance_bin, signature) %>%
    summarize(
      mean_score = mean(mean_score, na.rm = TRUE),
      se_score = sd(mean_score, na.rm = TRUE) / sqrt(n()),
      n_samples = n(),
      .groups = "drop"
    )

  p2 <- ggplot(condition_summary, aes(x = distance_bin, y = mean_score, color = condition, group = condition)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score), width = 0.2) +
    facet_wrap(~signature, scales = "free_y") +
    labs(title = "Signature Scores by Distance (Averaged by Condition)",
         x = "Distance (μm)", y = "Mean Score") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_figure(file.path(OUT_FIGURES_SPATIAL, "distance_signatures_by_condition.png"), p2, width = 12, height = 8)

  # Inflammation gradient (combined Implant_Up + DAM)
  if (all(c("Implant_Up", "DAM") %in% distance_summary$signature)) {
    inflammation_gradient <- distance_summary %>%
      filter(signature %in% c("Implant_Up", "DAM")) %>%
      group_by(sample, condition, distance_bin) %>%
      summarize(inflammation_score = mean(mean_score, na.rm = TRUE), .groups = "drop")

    p3 <- ggplot(inflammation_gradient, aes(x = distance_bin, y = inflammation_score, fill = condition)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.1, alpha = 0.5) +
      labs(title = "Inflammation Gradient from Electrode",
           subtitle = "Combined Implant_Up + DAM scores",
           x = "Distance (μm)", y = "Inflammation Score") +
      theme_publication()

    save_figure(file.path(OUT_FIGURES_SPATIAL, "distance_inflammation_gradient.png"), p3, width = 10, height = 6)
  }

  # Neuronal preservation by distance
  if ("Neuronal" %in% distance_summary$signature) {
    neuronal_gradient <- distance_summary %>%
      filter(signature == "Neuronal")

    p4 <- ggplot(neuronal_gradient, aes(x = distance_bin, y = mean_score, fill = condition)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.1, alpha = 0.5) +
      labs(title = "Neuronal Activity by Distance from Electrode",
           x = "Distance (μm)", y = "Neuronal Score") +
      theme_publication()

    save_figure(file.path(OUT_FIGURES_SPATIAL, "distance_neuronal_gradient.png"), p4, width = 10, height = 6)
  }
}

# Detailed spot-level analysis
if (length(all_distance_data) > 0) {
  spot_df <- bind_rows(all_distance_data)
  write.csv(spot_df, file.path(OUT_TABLES_SPATIAL, "distance_spot_data.csv"), row.names = FALSE)

  # Correlation: distance vs signature scores
  if ("Implant_Up" %in% colnames(spot_df) && "distance_um" %in% colnames(spot_df)) {
    cor_inflam <- cor.test(spot_df$distance_um, spot_df$Implant_Up, method = "spearman")
    cat(sprintf("\nDistance-Inflammation correlation: rho = %.3f (p = %.2e)\n",
                cor_inflam$estimate, cor_inflam$p.value))
  }

  if ("Neuronal" %in% colnames(spot_df) && "distance_um" %in% colnames(spot_df)) {
    cor_neur <- cor.test(spot_df$distance_um, spot_df$Neuronal, method = "spearman")
    cat(sprintf("Distance-Neuronal correlation: rho = %.3f (p = %.2e)\n",
                cor_neur$estimate, cor_neur$p.value))
  }
}

cat(sprintf("\nSaved distance analysis to: %s/\n", OUT_TABLES_SPATIAL))
