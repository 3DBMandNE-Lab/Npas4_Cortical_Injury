# config.R Reference

All analysis scripts `source("R/config.R")` as their first step. Scripts are organized in stage-based subdirectories under `R/` (e.g., `R/01_differential_expression/`, `R/04_spatial/`) but are always run from the project root, so the `source("R/config.R")` path remains unchanged. This file defines paths, thresholds, color palettes, a publication-quality ggplot2 theme, and a helper function for saving figures.

---

## Data Paths

| Variable | Default Value | Description |
|----------|--------------|-------------|
| `DATA_INTERNAL` | `"data/Internal"` | Root for internal datasets (Silicon RNA-seq, Polyimide Microarray) |
| `DATA_EXTERNAL` | `"data/external"` | Root for external datasets (Visium, human DBS) |
| `DATA_PROCESSED` | `"data/processed"` | Cached intermediate files (Seurat objects, metadata) |
| `SNRNASEQ_PATH` | `data/external/snrnaseq/seu_merged_annotated.RDS` | Path to the 19 GB snRNA-seq Seurat object (81,834 cells) |

### Overriding SNRNASEQ_PATH

The snRNA-seq Seurat object is large and may live on a different volume. Override it with an environment variable before running any 06* script:

```r
Sys.setenv(SNRNASEQ_PATH = "/your/path/seu_merged_annotated.RDS")
```

Or set it in your `.Renviron` file:

```
SNRNASEQ_PATH=/path/to/your/seu_merged_annotated.RDS
```

The config reads this with `Sys.getenv("SNRNASEQ_PATH", unset = <default>)`.

---

## Output Paths

| Variable | Path | Content |
|----------|------|---------|
| `OUT_TABLES` | `output/tables` | Root for all CSV results |
| `OUT_TABLES_DEG` | `output/tables/deg` | DESeq2 and limma results |
| `OUT_TABLES_COMPARISON` | `output/tables/comparison` | Cross-platform validation, validated genes |
| `OUT_TABLES_SPATIAL` | `output/tables/spatial` | Visium statistics, kill zone, correlations |
| `OUT_TABLES_SNRNASEQ` | `output/tables/snrnaseq` | snRNA-seq DE, signatures, microglia, neurons |
| `OUT_TABLES_HUMAN` | `output/tables/human` | Human DBS validation |
| `OUT_FIGURES` | `output/figures` | Root for all figures |
| `OUT_FIGURES_DEG` | `output/figures/deg` | Volcano plots, MA plots |
| `OUT_FIGURES_COMPARISON` | `output/figures/comparison` | Cross-platform scatter plots |
| `OUT_FIGURES_SPATIAL` | `output/figures/spatial` | Spatial heatmaps, distance plots |
| `OUT_FIGURES_SNRNASEQ` | `output/figures/snrnaseq` | Cell-type plots, violins |
| `OUT_FIGURES_HUMAN` | `output/figures/human` | Human validation plots |
| `OUT_FIGURES_MANUSCRIPT` | `output/figures/manuscript` | Publication-ready figures |

All directories are created automatically on `source("R/config.R")`.

---

## Analysis Thresholds

| Variable | Value | Usage |
|----------|-------|-------|
| `FDR_THRESH` | `0.05` | FDR cutoff for DE significance |
| `LFC_THRESH` | `0.5` | log2 fold change cutoff (used in pseudobulk DE) |

---

## Color Palettes

### Primary Categorical (4 categories)

`COL_CAT4`: `#E41A1C`, `#4DAF4A`, `#984EA3`, `#377EB8`

### Extended Categorical (8 categories)

`COL_CAT8`: `#E41A1C`, `#377EB8`, `#4DAF4A`, `#984EA3`, `#FF7F00`, `#FFFF33`, `#A65628`, `#F781BF`

### Sequential (low to high)

`COL_SEQ`: `#DEEBF7` --> `#9ECAE1` --> `#3182BD`

### Diverging

| Variable | Hex | Meaning |
|----------|-----|---------|
| `COL_DIV_NEG` | `#D73027` | Negative / downregulated |
| `COL_DIV_NEU` | `#FFFFBF` | Neutral / zero |
| `COL_DIV_POS` | `#1A9850` | Positive / upregulated |

### Utility Colors

| Variable | Hex | Usage |
|----------|-----|-------|
| `COL_SIG` | `#4393C3` | Significant points (steel blue) |
| `COL_NS` | `#BDBDBD` | Non-significant points (grey) |
| `COL_REF` | `#757575` | Reference lines (dark grey) |
| `COL_HIGHLIGHT` | `#FF7F00` | Highlighted elements (orange) |

### Direction Colors

| Variable | Hex | Usage |
|----------|-----|-------|
| `COL_UP` | `#E41A1C` | Upregulated / Implant |
| `COL_DOWN` | `#377EB8` | Downregulated / Stab |

### Condition Colors

`COL_CONDITION`: named vector

| Key | Hex | Color |
|-----|-----|-------|
| `"Control"` | `#4DAF4A` | Green |
| `"Implant"` | `#E41A1C` | Red |
| `"Stab"` | `#377EB8` | Blue |

Also available as individual variables: `COL_CONTROL`, `COL_IMPLANT`, `COL_STAB`.

### Spatial Signature Colors

`COL_SPATIAL`: named vector

| Key | Hex |
|-----|-----|
| `"Implant_Up"` | `#E41A1C` |
| `"DAM"` | `#984EA3` |
| `"Neuronal"` | `#4DAF4A` |
| `"Complement"` | `#377EB8` |

---

## Publication Theme

`theme_publication(base_size = 11, base_family = "")` returns a ggplot2 theme built on `theme_minimal` with:

- Black panel border (0.8pt)
- Subtle grey92 major gridlines, no minor gridlines
- White panel background
- Black axis ticks and text
- Bold title, grey subtitle, right-aligned caption
- White legend background
- Grey95 strip backgrounds for facets
- 10pt margins

---

## Figure Saving Helper

```r
save_figure(filename, plot, width = 7, height = 6, dpi = 300)
```

Saves both PNG (at specified dpi) and PDF (vector) versions. The PDF filename is derived by replacing `.png` with `.pdf`. Both are saved with a white background.
