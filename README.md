# A Conserved Multi-Lineage Transcriptional Architecture at the Neural Electrode Interface

Companion analysis code for [Töykkälä et al., 2026].

## Overview

Analysis pipeline for integrating bulk transcriptomics from two mechanically distinct electrode platforms (silicon and polyimide), single-nucleus RNA-seq, and Visium spatial transcriptomics to characterize the cellular programs at the neural electrode interface.

See `docs/` for per-script documentation.

## Quick Start

```bash
# 1. Install R dependencies
Rscript R/install_packages.R

# 2. Set path to snRNA-seq Seurat object (19 GB, not included)
export SNRNASEQ_PATH="/path/to/seu_merged_annotated.RDS"

# 3. Run analyses in order (see Pipeline below)
Rscript R/01_differential_expression/01_silicon_deseq2.R
```

## Data Availability

| Dataset | Accession | Description |
|---------|-----------|-------------|
| Silicon Probe | [GEO accession TBD] | 48 samples |
| Polyimide Probe | [GEO accession TBD] | 58 samples  |
| snRNA-seq | [GEO accession TBD] | 81,834 nuclei, 28 samples |
| Visium spatial | [GSE202425] | 6 tissue sections |

The snRNA-seq Seurat object is provided on zenodo [Zenodo accession TBD]. Set the `SNRNASEQ_PATH` environment variable to its location.

## Pipeline

All scripts run from the project root and source `R/config.R` for paths, thresholds, and palettes.

### Stage 1: Differential Expression
```
R/01_differential_expression/01_silicon_deseq2.R            # DESeq2 on silicon RNA-seq (Implant vs Control)
R/01_differential_expression/02_polyimide_limma.R           # limma on polyimide microarray (3-way comparison)
R/01_differential_expression/02b_timecourse_analysis.R      # Silicon temporal dynamics (1d-6wk)
R/01_differential_expression/02c_polyimide_timecourse.R     # Polyimide temporal dynamics
R/01_differential_expression/02d_concordant_timecourse.R    # Cross-platform temporal concordance
```

### Stage 2: Cross-Platform Concordance
```
R/02_cross_platform/03_cross_platform.R            # Validate polyimide signature in silicon
R/02_cross_platform/03a_cross_platform_stats.R     # Direction concordance, rank correlation
```

### Stage 3: Enrichment & TF Activity
```
R/03_enrichment/04b_gsea_enrichment.R          # GSEA pathway enrichment (KEGG, GO)
R/03_enrichment/04c_cell_deconvolution.R       # Cell-type deconvolution from bulk
R/03_enrichment/04d_tf_activity.R              # TF activity inference (dorothea/viper)
```

### Stage 4: Spatial Transcriptomics
```
R/04_spatial/05a_load_spatial.R             # Load Visium data (6 sections)
R/04_spatial/05b_score_spatial.R            # Score with concordant signatures
R/04_spatial/05c_spatial_stats.R            # Statistical analysis
R/04_spatial/05d_spatial_dimensions.R       # Injury zone dimensions
R/04_spatial/05e_spatial_morans_i.R         # Moran's I spatial autocorrelation
R/04_spatial/05e_spatial_killzone_summary.R # Kill zone summary statistics
R/04_spatial/05f_spatial_killzone.R         # Kill zone quantification
R/04_spatial/05g_spatial_correlations.R     # Pathway spatial correlations
R/04_spatial/05h_spatial_exclusion.R        # Domain exclusion analysis
R/04_spatial/05k_spatial_distance.R         # Distance-based analysis
R/04_spatial/05l_spatial_convergence.R      # SPP1/Complement/NPAS4 co-occurrence
R/04_spatial/05m_spatial_subtype_validation.R # Neuronal subtype spatial validation
```

### Stage 5: snRNA-seq Analysis (requires Seurat object)
```
R/05_snrnaseq/06a_load_snrnaseq.R            # Load dataset
R/05_snrnaseq/06a_celltype_attribution.R     # Cell-type attribution of concordant genes
R/05_snrnaseq/06b_score_snrnaseq.R           # Score cells with concordant signatures
R/05_snrnaseq/06b_stab_comparison.R          # Implant vs Stab wound comparison
R/05_snrnaseq/06c_snrnaseq_stats.R           # Statistical analysis
R/05_snrnaseq/06e_snrnaseq_pseudobulk.R      # Pseudobulk DE analysis
R/05_snrnaseq/06f_snrnaseq_signatures.R      # Signature scoring
R/05_snrnaseq/06g_microglial_states.R        # Microglial state analysis
R/05_snrnaseq/06h_spp1_microglia_analysis.R  # SPP1+ microglia characterization
R/05_snrnaseq/06i_complement_microglia_analysis.R  # Complement-high microglia
R/05_snrnaseq/06j_microglia_populations_summary.R  # Population overlap (Fisher's test)
R/05_snrnaseq/06k_spatial_spp1_complement.R  # Spatial colocalization of populations
R/05_snrnaseq/06l_microglia_temporal_pseudobulk.R  # Temporal pseudobulk
R/05_snrnaseq/06l_spp1_implant_vs_stab.R     # SPP1 implant-specificity test
R/05_snrnaseq/06m_spp1_proportion_test.R     # SPP1+ enrichment test
R/05_snrnaseq/06n_neuronal_pseudobulk_de.R   # Neuronal pseudobulk DE
R/05_snrnaseq/06n_neuronal_subtype_de.R      # Neuronal subtype DE (excitatory, PV+, CCK+)
R/05_snrnaseq/06o_oligo_subtype_de.R         # Oligodendrocyte subtype DE
R/05_snrnaseq/06o_oligodendrocyte_pseudobulk_de.R  # Oligodendrocyte pseudobulk DE
```

### Stage 6: Astrocyte Analysis
```
R/06_astrocyte/08_astrocyte_analysis.R        # Astrocyte state analysis
R/06_astrocyte/08b_astrocyte_spatial.R        # Astrocyte spatial patterns
```

### Stage 7: Robustness & Sensitivity
```
R/07_robustness/09a_spatial_density_permutation.R   # Density-controlled colocalization test
R/07_robustness/09b_killzone_dimensions_robust.R    # Bootstrap CIs for injury zone
R/07_robustness/09c_population_bootstrap.R          # Bootstrap for SPP1+ population
R/07_robustness/10a_concordance_sensitivity.R       # Concordance threshold sensitivity
R/07_robustness/10b_bootstrap_ci_check.R           # Bootstrap CI method validation
```

## Directory Structure

```
R/                              # Analysis scripts + config + installer
├── config.R                    # Shared paths, thresholds, palettes
├── install_packages.R          # Dependency installer
├── 01_differential_expression/ # DE analyses (01, 02, 02b-02d)
├── 02_cross_platform/          # Cross-platform concordance (03, 03a)
├── 03_enrichment/              # GSEA, deconvolution, TF activity (04b-04d)
├── 04_spatial/                 # Visium spatial analysis (05a-05m)
├── 05_snrnaseq/                # snRNA-seq analysis (06a-06o)
├── 06_astrocyte/               # Astrocyte analysis (08, 08b)
└── 07_robustness/              # Robustness & sensitivity (09a-10b)
docs/                           # Per-script documentation
data/                           # Input data (not tracked)
output/                         # Results (auto-created by scripts)
```

## Configuration

`R/config.R` defines all paths (relative to project root), analysis thresholds, color palettes, and a publication theme. See `docs/config.md` for details.

## Requirements

- R >= 4.2.0
- 13 Bioconductor + 14 CRAN packages (see `R/install_packages.R`)

## Citation

See `CITATION.cff` for citation metadata. If you use this code, please cite:

> [Töykkälä et al.] "A Conserved Multi-Lineage Transcriptional Architecture at the Neural Electrode Interface." , 2026.

## Corresponding Author

**Kevin Joseph**,

Laboratory for Neuroengineering, 

Department of Neurosurgery, 

Medical Center - University of Freiburg 

(ORCID: [0000-0001-6317-8736])

## License

This project is licensed under the MIT License — see `LICENSE` for details.
