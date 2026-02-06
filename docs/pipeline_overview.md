# Pipeline Overview

This document describes the full analysis pipeline for the neural electrode biocompatibility study comparing silicon and polyimide implants. Scripts are organized in numbered stages and should be executed in order, respecting the dependency structure below.

## Prerequisites

1. Run `R/install_packages.R` to install all CRAN and Bioconductor dependencies.
2. All scripts are organized in stage-based subdirectories under `R/` but still `source("R/config.R")` from the project root for paths, thresholds, color palettes, and the publication theme. Always run scripts from the project root directory.
3. The snRNA-seq Seurat object (19 GB) must be accessible at the path configured via `SNRNASEQ_PATH`.

---

## Pipeline Stages

### Stage 1: Differential Expression (01-02)

Raw data in, DE results out. No cross-script dependencies within this stage.

- **01** Silicon DESeq2 (RNA-seq counts --> DEGs)
- **02** Polyimide limma (CEL files --> DEGs + three-way comparison)
- **02b** Silicon timecourse (RNA-seq counts + metadata --> timepoint DEGs)
- **02c** Polyimide timecourse (CEL files --> timepoint DEGs)

### Stage 2: Cross-Platform Concordance (03)

Requires Stage 1 outputs.

- **03** Cross-platform validation (implant_specific + silicon DEGs --> validated_genes)
- **03a** Full concordance statistics (all DEGs + implant_specific --> concordance_statistics)

### Stage 2b: Concordant Timecourse (02d)

Requires Stages 1 and 2.

- **02d** Concordant gene timecourse (validated_genes + timecourse DEGs --> temporal trajectories)

### Stage 3: Enrichment and Deconvolution (04)

Requires Stage 1 and optionally Stage 2 outputs.

- **04b** GSEA enrichment (DEGs --> pathway NES, pathway concordance)
- **04c** Cell deconvolution (DEGs + snRNA-seq reference --> cell-type attribution)
- **04d** TF activity (DEGs --> TF activity scores, TF concordance)

### Stage 4: Spatial Transcriptomics (05)

Sequential dependency chain within this stage.

- **05a** Load Visium data (raw 10X --> cached Seurat list)
- **05b** Score spatial (cached Seurat + implant_specific --> scored Seurat)
- **05c-05h** Spatial statistics, dimensions, kill zone, correlations, exclusion (scored Seurat --> tables)
- **05e** Moran's I spatial autocorrelation (scored Seurat --> statistics)
- **05k** Distance-dependent analysis (scored Seurat --> distance gradients)
- **05l** Spatial convergence (raw Visium or cached Seurat --> convergence zones)
- **05m** Spatial subtype validation (cached Seurat --> subtype correlation with inflammation)

### Stage 5: snRNA-seq Analysis (06)

Requires Stage 2 outputs (validated_genes) and the external snRNA-seq Seurat object.

- **06a_load** Load and inspect snRNA-seq dataset
- **06a_celltype** Cell-type attribution of concordant genes
- **06b_score** Add module scores to snRNA-seq
- **06b_stab** Stab wound comparison + NPAS4 gradient
- **06c** Summary statistics
- **06e** Pseudobulk DE (all cells)
- **06f** Neuronal signature analysis
- **06g** Microglial state analysis
- **06h** SPP1+ microglia deep characterization
- **06i** Complement+ microglia deep characterization
- **06j** Microglial population summary (all populations)
- **06k** Spatial SPP1 vs Complement localization
- **06l_temporal** Microglia temporal dynamics (pseudobulk)
- **06l_spp1** SPP1 implant vs stab DE test
- **06m** SPP1+ proportion test
- **06n_pseudobulk** Neuronal pseudobulk DE
- **06n_subtype** Neuronal subtype-specific DE
- **06o_subtype** Oligodendrocyte subtype-specific DE
- **06o_pseudobulk** Oligodendrocyte pseudobulk DE

### Stage 6: Astrocyte Analysis (08)

Requires snRNA-seq data and spatial data.

- **08** Astrocyte snRNA-seq analysis (states, DE, complement production)
- **08b** Astrocyte spatial analysis (two-cell complement model)

### Stage 7: Robustness Analyses (09)

Requires Stages 4 and 5 outputs.

- **09a** Density-controlled spatial permutation test
- **09b** Robust kill zone dimensions with bootstrap CIs
- **09c** Population bootstrap CIs (microglial overlap statistics)

### Stage 8: Sensitivity Analyses (10)

Requires Stages 2 and 5 outputs.

- **10a** Concordance sensitivity at relaxed thresholds
- **10b** Bootstrap CI method verification

---

## Dependency Diagram

```
DATA INPUTS
============
Silicon RNA-seq counts ─────────┐
Polyimide CEL files ────────────┤
Silicon metadata (Excel) ───────┤
Visium 10X matrices ────────────┤
snRNA-seq Seurat (19 GB) ──────┤
                                │
STAGE 1: DE (01_differential_expression/)
=========                       │
01_silicon_deseq2 ──────────────┤──► deg/silicon_deseq2_results.csv
02_polyimide_limma ─────────────┤──► deg/polyimide_limma_results.csv
                                │    comparison/implant_specific.csv
02b_timecourse_analysis ────────┤──► timecourse/de_*.csv
02c_polyimide_timecourse ───────┤──► timecourse/polyimide_*.csv
                                │
STAGE 2: CONCORDANCE (02_cross_platform/)
====================            │
03_cross_platform ──────────────┤──► comparison/validated_genes.csv
    (needs: 01, 02)             │    comparison/cross_platform_validation.csv
03a_cross_platform_stats ───────┤──► comparison/concordance_statistics.csv
    (needs: 01, 02)             │    comparison/implant_signature_validation_stats.csv
                                │
02d_concordant_timecourse ──────┤──► timecourse/concordant_signature_scores.csv
    (in 01_differential_expression/)    timecourse/concordant_genes_timecourse.csv
    (needs: 02b, 02c, 03)      │
                                │
STAGE 3: ENRICHMENT (03_enrichment/)
===================             │
04b_gsea_enrichment ────────────┤──► enrichment/gsea_*.csv, pathway_concordance.csv
    (needs: 01, 02)             │
04c_cell_deconvolution ─────────┤──► deconvolution/concordant_celltype_attribution.csv
    (needs: 03, snRNA-seq)      │
04d_tf_activity ────────────────┤──► enrichment/tf_activity_*.csv, tf_concordance.csv
    (needs: 01, 02)             │
                                │
STAGE 4: SPATIAL (04_spatial/)   │
================                │
05a_load_spatial ───────────────┤──► data/processed/spatial_seurat_list.RDS
    (needs: Visium data)        │
05b_score_spatial ──────────────┤──► data/processed/spatial_scored.RDS
    (needs: 05a, 02)            │
05c_spatial_stats ──────────────┤──► spatial/spatial_stats.csv
05d_spatial_dimensions ─────────┤──► spatial/killzone_dimensions.csv
05e_spatial_killzone_summary ───┤──► (figures only)
05e_spatial_morans_i ───────────┤──► spatial/morans_i_statistics.csv
05f_spatial_killzone ───────────┤──► spatial/killzone_statistics.csv
05g_spatial_correlations ───────┤──► spatial/pathway_correlations.csv
05h_spatial_exclusion ──────────┤──► spatial/domain_overlap_matrix.csv
    (all 05c-05h need: 05b)     │
05k_spatial_distance ───────────┤──► spatial/distance_signature_summary.csv
05l_spatial_convergence ────────┤──► spatial/convergence/convergence_statistics.csv
05m_spatial_subtype_validation ─┤──► spatial/subtype_validation/*.csv
    (needs: 05a)                │
                                │
STAGE 5: snRNA-seq (05_snrnaseq/)
==================              │
06a_load_snrnaseq ──────────────┤──► data/processed/snrnaseq_meta.csv
06a_celltype_attribution ───────┤──► snrnaseq/concordant_celltype_attribution.csv
    (needs: 03, snRNA-seq)      │
06b_score_snrnaseq ─────────────┤──► (module scores in memory)
06b_stab_comparison ────────────┤──► stab_comparison/concordant_gene_classification.csv
    (needs: 03, 02b, snRNA-seq) │
06c_snrnaseq_stats ─────────────┤──► snrnaseq/snrnaseq_summary.csv
06e_snrnaseq_pseudobulk ────────┤──► snrnaseq/pseudobulk_DE_*.csv
06f_snrnaseq_signatures ────────┤──► snrnaseq/signature_scores_summary.csv
06g_microglial_states ──────────┤──► snrnaseq/microglia_state_distribution.csv
06h_spp1_microglia_analysis ────┤──► snrnaseq/spp1_microglia/*.csv
06i_complement_microglia ───────┤──► snrnaseq/complement_microglia/*.csv
06j_microglia_populations ──────┤──► snrnaseq/microglia_populations/*.csv
06k_spatial_spp1_complement ────┤──► spatial/spp1_complement/*.csv
    (needs: 05a or cache)       │
06l_microglia_temporal ─────────┤──► snrnaseq/microglia_temporal/*.csv
06l_spp1_implant_vs_stab ───────┤──► snrnaseq/spp1_specificity/implant_vs_stab_de.csv
06m_spp1_proportion_test ───────┤──► snrnaseq/spp1_specificity/spp1_proportion_tests.csv
06n_neuronal_pseudobulk_de ─────┤──► snrnaseq/neurons/de_implant_vs_control.csv
06n_neuronal_subtype_de ────────┤──► snrnaseq/neurons/de_*_implant_vs_ctrl.csv
06o_oligo_subtype_de ───────────┤──► snrnaseq/oligodendrocytes/de_*_implant_vs_ctrl.csv
06o_oligodendrocyte_pseudobulk ─┤──► snrnaseq/oligodendrocytes/de_implant_vs_control.csv
    (all 06* need: snRNA-seq)   │
                                │
STAGE 6: ASTROCYTES (06_astrocyte/)
===================             │
08_astrocyte_analysis ──────────┤──► snrnaseq/astrocytes/*.csv
    (needs: snRNA-seq)          │
08b_astrocyte_spatial ──────────┤──► spatial/astrocyte/*.csv
    (needs: Visium data)        │
                                │
STAGE 7: ROBUSTNESS (07_robustness/)
===================             │
09a_spatial_density_permutation ┤──► spatial/permutation/density_permutation_results.csv
    (needs: 05b or 05a)        │
09b_killzone_dimensions_robust ─┤──► spatial/dimensions/killzone_dimensions_robust.csv
    (needs: 05b or 05a)        │
09c_population_bootstrap ───────┤──► snrnaseq/population_bootstrap/*.csv
    (needs: snRNA-seq)          │
                                │
STAGE 8: SENSITIVITY (07_robustness/)
====================            │
10a_concordance_sensitivity ────┤──► comparison/concordance_sensitivity.csv
    (needs: 03a)                │
10b_bootstrap_ci_check ─────────┘──► snrnaseq/bootstrap_ci_comparison.csv
    (needs: 06f)
```

---

## Quick Start

```r
# 1. Install dependencies
source("R/install_packages.R")

# 2. Run DE analyses
source("R/01_differential_expression/01_silicon_deseq2.R")
source("R/01_differential_expression/02_polyimide_limma.R")

# 3. Cross-platform validation
source("R/02_cross_platform/03_cross_platform.R")
source("R/02_cross_platform/03a_cross_platform_stats.R")

# 4. Enrichment (can run in parallel)
source("R/03_enrichment/04b_gsea_enrichment.R")
source("R/03_enrichment/04d_tf_activity.R")

# 5. Spatial analysis (sequential)
source("R/04_spatial/05a_load_spatial.R")
source("R/04_spatial/05b_score_spatial.R")
source("R/04_spatial/05c_spatial_stats.R")
# ... etc.

# 6. snRNA-seq analyses (each script loads its own data)
source("R/05_snrnaseq/06a_celltype_attribution.R")
# ... etc.
```
