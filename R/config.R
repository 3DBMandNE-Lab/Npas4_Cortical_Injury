# config.R - Paths, constants, and publication theme

# Data paths
DATA_INTERNAL <- "data/Internal"
DATA_EXTERNAL <- "data/external"
DATA_PROCESSED <- "data/processed"

# External snRNA-seq Seurat object (19 GB)
# Override with: Sys.setenv(SNRNASEQ_PATH = "/your/path/seu_merged_annotated.RDS")
SNRNASEQ_PATH <- Sys.getenv(
  "SNRNASEQ_PATH",
  unset = file.path(DATA_EXTERNAL, "snrnaseq", "seu_merged_annotated.RDS")
)

# Output paths - organized by analysis type
OUT_TABLES <- "output/tables"
OUT_TABLES_DEG <- file.path(OUT_TABLES, "deg")
OUT_TABLES_COMPARISON <- file.path(OUT_TABLES, "comparison")
OUT_TABLES_SPATIAL <- file.path(OUT_TABLES, "spatial")
OUT_TABLES_SNRNASEQ <- file.path(OUT_TABLES, "snrnaseq")
OUT_TABLES_HUMAN <- file.path(OUT_TABLES, "human")

OUT_FIGURES <- "output/figures"
OUT_FIGURES_DEG <- file.path(OUT_FIGURES, "deg")
OUT_FIGURES_COMPARISON <- file.path(OUT_FIGURES, "comparison")
OUT_FIGURES_SPATIAL <- file.path(OUT_FIGURES, "spatial")
OUT_FIGURES_SNRNASEQ <- file.path(OUT_FIGURES, "snrnaseq")
OUT_FIGURES_HUMAN <- file.path(OUT_FIGURES, "human")
OUT_FIGURES_MANUSCRIPT <- file.path(OUT_FIGURES, "manuscript")

# Create all directories
for (d in c(DATA_PROCESSED,
            OUT_TABLES_DEG, OUT_TABLES_COMPARISON, OUT_TABLES_SPATIAL,
            OUT_TABLES_SNRNASEQ, OUT_TABLES_HUMAN,
            OUT_FIGURES_DEG, OUT_FIGURES_COMPARISON, OUT_FIGURES_SPATIAL,
            OUT_FIGURES_SNRNASEQ, OUT_FIGURES_HUMAN, OUT_FIGURES_MANUSCRIPT)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Thresholds
FDR_THRESH <- 0.05
LFC_THRESH <- 0.5

# ============================================================
# COLOR PALETTES (from style guide)
# ============================================================

# Primary categorical (4 categories)
COL_CAT4 <- c("#E41A1C", "#4DAF4A", "#984EA3", "#377EB8")

# Extended categorical (6+ categories)
COL_CAT8 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
              "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

# Sequential (low to high)
COL_SEQ <- c("#DEEBF7", "#9ECAE1", "#3182BD")

# Diverging (negative / neutral / positive)
COL_DIV_NEG <- "#D73027"
COL_DIV_NEU <- "#FFFFBF"
COL_DIV_POS <- "#1A9850"

# Utility colors
COL_SIG <- "#4393C3"      # Significant (steel blue)
COL_NS <- "#BDBDBD"       # Non-significant (grey)
COL_REF <- "#757575"      # Reference line (dark grey)
COL_HIGHLIGHT <- "#FF7F00" # Highlight (orange)

# Legacy colors for compatibility
COL_UP <- "#E41A1C"
COL_DOWN <- "#377EB8"

# Condition colors
COL_CONTROL <- "#4DAF4A"
COL_IMPLANT <- "#E41A1C"
COL_STAB <- "#377EB8"
COL_CONDITION <- c("Control" = COL_CONTROL, "Implant" = COL_IMPLANT, "Stab" = COL_STAB)

# Spatial signature colors
COL_SPATIAL <- c(
  "Implant_Up" = "#E41A1C",
  "DAM" = "#984EA3",
  "Neuronal" = "#4DAF4A",
  "Complement" = "#377EB8"
)

# ============================================================
# PUBLICATION THEME (from style guide)
# ============================================================

theme_publication <- function(base_size = 11, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # Panel
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),

      # Axes
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.4),
      axis.text = ggplot2::element_text(color = "black", size = base_size - 1),
      axis.title = ggplot2::element_text(color = "black", size = base_size),

      # Text
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 3, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size, color = "grey30", hjust = 0),
      plot.caption = ggplot2::element_text(size = base_size - 2, color = "grey50", hjust = 1),

      # Legend
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      legend.title = ggplot2::element_text(face = "bold", size = base_size),
      legend.text = ggplot2::element_text(size = base_size - 1),

      # Margins
      plot.margin = ggplot2::margin(10, 12, 10, 10),

      # Strip (for facets)
      strip.background = ggplot2::element_rect(fill = "grey95", color = "black", linewidth = 0.4),
      strip.text = ggplot2::element_text(face = "bold", size = base_size)
    )
}

# Helper for saving figures
save_figure <- function(filename, plot, width = 7, height = 6, dpi = 300) {
  ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi, bg = "white")
  # Also save PDF version
  pdf_name <- sub("\\.png$", ".pdf", filename)
  ggplot2::ggsave(pdf_name, plot, width = width, height = height, bg = "white")
}
