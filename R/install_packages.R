# install_packages.R — Install all dependencies for the analysis pipeline
#
# Usage: Rscript R/install_packages.R

message("Installing CRAN packages...")
cran_packages <- c(
  "tidyverse", "readr", "readxl", "dplyr", "tidyr", "stringr",
  "ggplot2", "ggrepel", "patchwork", "scales", "viridis", "pheatmap",
  "Seurat", "SeuratObject",
  "openxlsx", "circlize", "here"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("  Installing: ", pkg)
    install.packages(pkg, quiet = TRUE)
  } else {
    message("  OK: ", pkg)
  }
}

message("\nInstalling Bioconductor packages...")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", quiet = TRUE)
}

bioc_packages <- c(
  "DESeq2", "limma", "edgeR",
  "clusterProfiler", "org.Rn.eg.db", "enrichplot",
  "dorothea", "viper",
  "oligo", "affycoretools",
  "pd.clariom.s.rat", "clariomsrathttranscriptcluster.db",
  "ComplexHeatmap"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("  Installing: ", pkg)
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    message("  OK: ", pkg)
  }
}

message("\nAll packages installed.")
message("Run sessionInfo() to check versions.")
