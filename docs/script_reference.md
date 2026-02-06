# Script Reference

Reference documentation for all analysis scripts in `R/`. Each entry documents the script's purpose, inputs, outputs, key parameters, and dependencies.

---

## Stage 01-02: Differential Expression (`R/01_differential_expression/`)

---

### 01_differential_expression/01_silicon_deseq2.R

**Purpose:** DESeq2 differential expression analysis of silicon (stiff) electrode RNA-seq data. Compares 6 implanted samples (1 week, 100um) vs 12 naive controls.

**Inputs:**
- `data/Internal/Silicon_RNAseq/RS1_near_far_deseq2_raw_counts.txt` (implant counts)
- `data/Internal/Silicon_RNAseq/Naivecontrols_deseq2_raw_counts.txt` (naive counts)

**Outputs:**
- `output/tables/deg/silicon_deseq2_results.csv`

**Key parameters:**
- FDR threshold: 0.05 (from config.R)
- Filter: genes with >= 10 counts in >= 2 samples
- Contrast: Implanted vs Control
- Implant samples: 6 (1 week, 100um); Control samples: 12 naive

**Dependencies:** None (first script in pipeline)

---

### 01_differential_expression/02_polyimide_limma.R

**Purpose:** limma differential expression analysis of polyimide (soft) electrode microarray data. Runs both a primary Implant vs Control comparison and a three-way comparison (Implant/Stab/Control) to identify implant-specific genes.

**Inputs:**
- `data/Internal/Polyimide_Microarray/*.CEL` (Clariom S Rat arrays)

**Outputs:**
- `output/tables/deg/polyimide_limma_results.csv` (Implant vs Control)
- `output/tables/comparison/implant_specific.csv` (genes sig in Implant vs Control but not Stab vs Control)
- `output/tables/comparison/shared_injury.csv` (genes sig in both)
- `output/tables/comparison/polyimide_three_way.csv` (full three-way results)

**Key parameters:**
- FDR threshold: 0.05
- Platform: Clariom S Rat (pd.clariom.s.rat)
- Normalization: RMA
- Gene aggregation: mean across probes per gene

**Dependencies:** None

---

### 01_differential_expression/02b_timecourse_analysis.R

**Purpose:** Timecourse DESeq2 analysis of silicon RNA-seq across multiple timepoints (24h, 1wk, 6wk) and distances (100um, 500um) vs naive control.

**Inputs:**
- `data/Internal/Silicon_RNAseq/RS1_near_far_deseq2_raw_counts.txt`
- `data/Internal/Silicon_RNAseq/Naivecontrols_deseq2_raw_counts.txt`
- `data/Internal/Silicon_RNAseq/Purcell RNAseq Sample Identification Key.xlsx`

**Outputs:**
- `output/tables/timecourse/de_24h_100um.csv`, `de_1wk_100um.csv`, etc. (6 pairwise comparisons)
- `output/tables/timecourse/deg_counts_by_condition.csv`
- `output/tables/timecourse/key_gene_expression.csv`
- `output/figures/timecourse/complement_trajectory.png`
- `output/figures/timecourse/neuronal_trajectory.png`
- `output/figures/timecourse/deg_heatmap.png`
- `output/figures/timecourse/killzone_trajectory.png`

**Key parameters:**
- 6 pairwise comparisons: each timepoint x distance vs Control
- Key genes tracked: Npas4, Arc, Fos, C1qa, C3, Spp1, Gfap, etc.

**Dependencies:** None

---

### 01_differential_expression/02c_polyimide_timecourse.R

**Purpose:** Timecourse limma analysis of polyimide microarray data across Week0, Week1, Week2, Week4, Week18 vs naive control.

**Inputs:**
- `data/Internal/Polyimide_Microarray/Trepanation_Implants/*.CEL`
- `data/Internal/Polyimide_Microarray/Naive_Controls/*.CEL`

**Outputs:**
- `output/tables/timecourse/polyimide_Week*_vs_Ctrl.csv` (5 comparisons)
- `output/tables/timecourse/polyimide_deg_counts.csv`
- `output/tables/timecourse/polyimide_key_gene_expression.csv`
- `output/figures/timecourse/polyimide_complement_trajectory.png`
- `output/figures/timecourse/polyimide_neuronal_trajectory.png`
- `output/figures/timecourse/polyimide_deg_heatmap.png`

**Key parameters:**
- 5 contrasts: each Week vs Control
- Platform: Clariom S Rat HT
- Normalization: RMA, max-mean probe selection per gene

**Dependencies:** Optionally reads `output/tables/timecourse/de_1wk_100um.csv` for cross-platform comparison (from 02b)

---

### 01_differential_expression/02d_concordant_timecourse.R

**Purpose:** Track the 50 validated concordant genes across both silicon and polyimide timecourse datasets. Calculates composite signature scores and percent-significant over time.

**Inputs:**
- `output/tables/comparison/validated_genes.csv` (from 03/03a)
- `output/tables/timecourse/de_*.csv` (from 02b)
- `output/tables/timecourse/polyimide_Week*.csv` (from 02c)

**Outputs:**
- `output/tables/timecourse/concordant_signature_scores.csv`
- `output/tables/timecourse/concordant_genes_timecourse.csv`
- `output/figures/timecourse/concordant_signature_trajectory.png`
- `output/figures/timecourse/concordant_up_trajectory.png`
- `output/figures/timecourse/concordant_pct_sig_trajectory.png`
- `output/figures/timecourse/concordant_key_genes_heatmap.png`
- `output/figures/timecourse/concordant_category_trajectory.png`
- `output/figures/timecourse/concordant_top_dynamic.png`

**Key parameters:**
- Signature score = mean log2FC of UP genes minus mean log2FC of DOWN genes
- Time alignment: silicon (24h, 1wk, 6wk) and polyimide (Week0-Week18)

**Dependencies:** 02b, 02c, 03 or 03a

---

## Stage 03: Cross-Platform Concordance (`R/02_cross_platform/`)

---

### 02_cross_platform/03_cross_platform.R

**Purpose:** Validate the polyimide implant-specific gene signature in the silicon dataset. Matches genes by name (case-insensitive), tests direction concordance, and identifies validated concordant genes.

**Inputs:**
- `output/tables/comparison/implant_specific.csv` (from 02)
- `output/tables/deg/silicon_deseq2_results.csv` (from 01)

**Outputs:**
- `output/tables/comparison/cross_platform_validation.csv`
- `output/tables/comparison/validated_genes.csv` (50 concordant genes)

**Key parameters:**
- Concordance criterion: significant in silicon (FDR < 0.05) AND same direction as polyimide
- Binomial test for direction concordance vs 50% chance

**Dependencies:** 01, 02

---

### 02_cross_platform/03a_cross_platform_stats.R

**Purpose:** Full cross-platform concordance statistics including genome-wide and implant-specific direction concordance, rank correlations, and validation metrics.

**Inputs:**
- `output/tables/deg/silicon_deseq2_results.csv` (from 01)
- `output/tables/deg/polyimide_limma_results.csv` (from 02)
- `output/tables/comparison/implant_specific.csv` (from 02)

**Outputs:**
- `output/tables/comparison/concordance_statistics.csv` (83.3% implant-assoc, 56.9% genome-wide)
- `output/tables/comparison/implant_signature_validation_stats.csv` (824 implant-specific, 557 matched, 50 validated)
- `output/tables/comparison/validated_genes.csv` (overwrites with same content)
- `output/tables/comparison/cross_platform_validation.csv`

**Key parameters:**
- Reports BOTH implant-associated (83.3%) and genome-wide (56.9%) concordance
- Uses Spearman rank correlation (valid across platforms with different dynamic ranges)

**Dependencies:** 01, 02

---

## Stage 04: Enrichment, Deconvolution, TF Activity (`R/03_enrichment/`)

---

### 03_enrichment/04b_gsea_enrichment.R

**Purpose:** Gene Set Enrichment Analysis (GSEA) using GO terms for both platforms. Calculates pathway-level concordance via NES correlation.

**Inputs:**
- `output/tables/deg/silicon_deseq2_results.csv` (from 01)
- `output/tables/deg/polyimide_limma_results.csv` (from 02)

**Outputs:**
- `output/tables/enrichment/gsea_silicon.csv`
- `output/tables/enrichment/gsea_polyimide.csv`
- `output/tables/enrichment/pathway_concordance.csv`
- `output/figures/enrichment/gsea_nes_correlation.png`
- `output/figures/enrichment/gsea_silicon_dotplot.png`
- `output/figures/enrichment/gsea_polyimide_dotplot.png`
- `output/figures/enrichment/gsea_top_pathways.png`

**Key parameters:**
- Ontology: GO ALL (BP + CC + MF)
- p-value cutoff: 0.25 (GSEA standard)
- Organism: org.Rn.eg.db (rat)

**Dependencies:** 01, 02

---

### 03_enrichment/04c_cell_deconvolution.R

**Purpose:** Cell-type deconvolution of concordant genes using the snRNA-seq reference. Determines which cell types express each concordant gene.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)
- `output/tables/comparison/validated_genes.csv` (from 03)

**Outputs:**
- `output/tables/deconvolution/celltype_markers.csv`
- `output/tables/deconvolution/concordant_celltype_attribution.csv`
- `output/figures/deconvolution/concordant_attribution_pie.png`
- `output/figures/deconvolution/concordant_neuronal_glial.png`

**Key parameters:**
- Marker detection: min.pct = 0.25, logfc.threshold = 0.5, top 50 per cell type
- Classification: >70% neuronal = "Neuronal", >70% glial = "Glial", else "Mixed"

**Dependencies:** 03, snRNA-seq data

---

### 03_enrichment/04d_tf_activity.R

**Purpose:** Transcription factor activity inference using DoRothEA regulons and VIPER. Computes TF activity concordance across platforms.

**Inputs:**
- `output/tables/deg/silicon_deseq2_results.csv` (from 01)
- `output/tables/deg/polyimide_limma_results.csv` (from 02)

**Outputs:**
- `output/tables/enrichment/tf_activity_silicon.csv`
- `output/tables/enrichment/tf_activity_polyimide.csv`
- `output/tables/enrichment/tf_concordance.csv`
- `output/figures/enrichment/tf_activity_concordance.png`
- `output/figures/enrichment/tf_top_activities.png`
- `output/figures/enrichment/tf_key_activities.png`

**Key parameters:**
- DoRothEA regulons: mouse (dorothea_mm), confidence A-C
- VIPER minimum regulon size: 5
- Key inflammatory TFs: Nfkb1, Rela, Stat1, Stat3, Irf1, Irf7, Jun, Fos
- Key neuronal TFs: Creb1, Mef2a, Mef2c, Neurod1, Rest, Nfia

**Dependencies:** 01, 02

---

## Stage 05: Spatial Transcriptomics (`R/04_spatial/`)

---

### 04_spatial/05a_load_spatial.R

**Purpose:** Load Visium 10X spatial transcriptomics data, normalize, add spatial coordinates, and cache as a list of Seurat objects.

**Inputs:**
- `data/external/stRNAseq/sample_sheet.csv`
- `data/external/stRNAseq/{sample_id}/filtered_feature_bc_matrix/`
- `data/external/stRNAseq/{sample_id}/spatial/tissue_positions.csv`

**Outputs:**
- `data/processed/spatial_seurat_list.RDS`

**Key parameters:**
- 6 Visium sections
- Normalization: Seurat NormalizeData (log-normalization)

**Dependencies:** None (raw data)

---

### 04_spatial/05b_score_spatial.R

**Purpose:** Add module scores (Implant_Up, DAM, Neuronal, Complement) to all spatial Seurat objects.

**Inputs:**
- `data/processed/spatial_seurat_list.RDS` (from 05a)
- `output/tables/comparison/implant_specific.csv` (from 02)

**Outputs:**
- `data/processed/spatial_scored.RDS`

**Key parameters:**
- Signatures: Implant_Up (implant-specific genes), DAM (Spp1, Tyrobp, Cd68, Aif1, Trem2, Apoe), Neuronal (Npas4, Arc, Fos, Egr1, Nr4a1, Bdnf), Complement (C1qa, C1qb, C1qc, C3)
- Scoring: Seurat AddModuleScore, nbin=10

**Dependencies:** 05a, 02

---

### 04_spatial/05c_spatial_stats.R

**Purpose:** Calculate basic spatial statistics including mean signature scores, kill zone percentage, and inter-signature correlations per sample.

**Inputs:**
- `data/processed/spatial_scored.RDS` (from 05b)

**Outputs:**
- `output/tables/spatial/spatial_stats.csv`

**Key parameters:**
- Kill zone: Implant_Up > Q75 AND Neuronal < Q25

**Dependencies:** 05b

---

### 04_spatial/05d_spatial_dimensions.R

**Purpose:** Quantify kill zone physical dimensions from Visium data using spot coordinates and 10X Visium specifications.

**Inputs:**
- `data/processed/spatial_scored.RDS` (from 05b)

**Outputs:**
- `output/tables/spatial/killzone_dimensions.csv`
- `output/figures/manuscript/fig_spatial_dimensions.png`

**Key parameters:**
- Visium spot diameter: 55 um
- Spot spacing: 100 um center-to-center
- Equivalent radius: sqrt(n_spots * 10000 / pi)
- Kill zone threshold: Implant_Up > Q75 AND Neuronal < Q25

**Dependencies:** 05b

---

### 04_spatial/05e_spatial_killzone_summary.R

**Purpose:** Summarize kill zone dimensions and compare Visium-measured radii to bulk RNA-seq sampling distances (100um vs 500um).

**Inputs:**
- `output/tables/spatial/killzone_dimensions.csv` (from 05d)

**Outputs:**
- `output/figures/manuscript/fig_killzone_vs_bulk_sampling.png`
- `output/figures/manuscript/fig_killzone_expansion.png`

**Key parameters:**
- Spot spacing: 100 um
- Reference distances: 100um (catastrophic zone), 500um (protected zone)

**Dependencies:** 05d

---

### 04_spatial/05e_spatial_morans_i.R

**Purpose:** Calculate Moran's I spatial autocorrelation statistic for each signature in each sample. Tests whether signature scores are spatially clustered.

**Inputs:**
- `data/processed/spatial_scored.RDS` (from 05b)

**Outputs:**
- `output/tables/spatial/morans_i_statistics.csv`

**Key parameters:**
- k = 6 nearest neighbors
- Manual implementation (no spdep dependency)
- Two-sided z-test for significance

**Dependencies:** 05b

---

### 04_spatial/05f_spatial_killzone.R

**Purpose:** Kill zone analysis quantifying spots with high inflammation and low neuronal activity, plus Implant_Up vs Neuronal correlation.

**Inputs:**
- `data/processed/spatial_scored.RDS` (from 05b)

**Outputs:**
- `output/tables/spatial/killzone_statistics.csv`

**Key parameters:**
- Kill zone: Implant_Up > Q75 AND Neuronal < Q25
- Pearson correlation between Implant_Up and Neuronal scores

**Dependencies:** 05b

---

### 04_spatial/05g_spatial_correlations.R

**Purpose:** Calculate all pairwise Pearson correlations between spatial signature scores across samples.

**Inputs:**
- `data/processed/spatial_scored.RDS` (from 05b)

**Outputs:**
- `output/tables/spatial/pathway_correlations.csv` (long format)
- `output/tables/spatial/pathway_correlations_wide.csv` (wide format)

**Key parameters:**
- Signatures: Implant_Up, DAM, Neuronal, Complement
- Pearson correlation with p-values

**Dependencies:** 05b

---

### 04_spatial/05h_spatial_exclusion.R

**Purpose:** Spatial domain overlap/exclusion analysis. Calculates Jaccard indices, overlap coefficients, fold enrichment, and Fisher's exact tests for all pairwise domain combinations.

**Inputs:**
- `data/processed/spatial_scored.RDS` (from 05b)

**Outputs:**
- `output/tables/spatial/domain_overlap_matrix.csv`
- `output/tables/spatial/domain_exclusion_summary.csv`

**Key parameters:**
- Domain threshold: top 25% for inflammatory signatures, bottom 25% for Neuronal
- Fold enrichment: observed overlap / expected overlap under independence

**Dependencies:** 05b

---

### 04_spatial/05k_spatial_distance.R

**Purpose:** Distance-dependent analysis of signature scores from estimated electrode center. Bins spots by distance and calculates mean scores per bin.

**Inputs:**
- `data/processed/spatial/{seurat_list.RDS}` (note: uses different path than other 05* scripts)

**Outputs:**
- `output/tables/spatial/distance_signature_summary.csv`
- `output/tables/spatial/distance_spot_data.csv`
- `output/figures/spatial/distance_signatures_by_sample.png`
- `output/figures/spatial/distance_signatures_by_condition.png`
- `output/figures/spatial/distance_inflammation_gradient.png`
- `output/figures/spatial/distance_neuronal_gradient.png`

**Key parameters:**
- Electrode center: centroid of top 5% inflammation spots
- Distance bins: 0-100, 100-200, 200-500, 500-1000, >1000 um
- Pixel-to-um conversion from spot spacing

**Dependencies:** 05a (or scored variant)

---

### 04_spatial/05l_spatial_convergence.R

**Purpose:** Spatial proximity analysis testing whether SPP1, Complement, and NPAS4-silenced domains converge in the same tissue zones. Classifies spots by domain membership and calculates enrichment statistics.

**Inputs:**
- `data/processed/spatial_seurat_list.RDS` (from 05a, or loads from raw)

**Outputs:**
- `output/tables/spatial/convergence/convergence_statistics.csv`
- `data/processed/spatial_convergence_annotated.RDS`
- `output/figures/manuscript/spatial_convergence/{sample}_convergence_maps.png`
- `output/figures/manuscript/spatial_convergence/{sample}_domain_distribution.png`
- `output/figures/manuscript/spatial_convergence/enrichment_heatmap.png`
- `output/figures/manuscript/spatial_convergence/convergence_by_condition.png`
- `output/figures/manuscript/spatial_convergence/complement_vs_npas4_scatter.png`

**Key parameters:**
- SPP1 genes: Spp1, Gpnmb, Lgals3, Fabp5, Igf1, Cd63, Lpl, Mmp12
- Complement genes: C1qa, C1qb, C1qc, C3
- Neuronal activity genes: Npas4, Arc, Fos, Egr1, Nr4a1, Bdnf
- Domain threshold: top/bottom 25%
- Convergence zone: spots meeting 2+ criteria
- Main display sample: Visium_7C (Chronic No Stim)

**Dependencies:** 05a

---

### 04_spatial/05m_spatial_subtype_validation.R

**Purpose:** Spatial validation of neuronal subtype vulnerability and oligodendrocyte findings. Tests whether excitatory markers show stronger spatial depletion near electrodes than inhibitory markers, and whether Qk is spatially depleted.

**Inputs:**
- `data/processed/spatial_seurat_list.RDS` (from 05a)

**Outputs:**
- `output/tables/spatial/subtype_validation/spatial_subtype_scores.csv`
- `output/tables/spatial/subtype_validation/spatial_subtype_correlations.csv`
- `output/tables/spatial/subtype_validation/spatial_correlation_summary.csv`
- `output/figures/manuscript/spatial_subtype_validation.png`

**Key parameters:**
- Excitatory signature: Camk2a, Grin2a, Slc17a7, Arc, Snap25, Syt1
- Inhibitory signature: Gad1, Gad2, Slc32a1, Pvalb, Sst
- Premyelinating: Qk (single gene)
- Paired Wilcoxon test comparing excitatory vs inhibitory inflammation correlations

**Dependencies:** 05a

---

## Stage 06: snRNA-seq Analysis (`R/05_snrnaseq/`)

---

### 05_snrnaseq/06a_load_snrnaseq.R

**Purpose:** Load and inspect the full snRNA-seq dataset. Saves metadata to disk for quick inspection without reloading the full object.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `data/processed/snrnaseq_meta.csv`

**Dependencies:** snRNA-seq data

---

### 05_snrnaseq/06a_celltype_attribution.R

**Purpose:** Cell-type attribution of the 50 concordant genes. Determines what fraction of each gene's expression comes from neurons, microglia, astrocytes, and oligodendrocytes.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)
- `output/tables/comparison/validated_genes.csv` (from 03)

**Outputs:**
- `output/tables/snrnaseq/concordant_celltype_attribution.csv`
- `output/tables/snrnaseq/concordant_expression_by_celltype.csv`
- `output/figures/manuscript/fig_celltype_attribution_bar.png`
- `output/figures/manuscript/fig_celltype_attribution_pie.png`
- `output/figures/manuscript/fig_celltype_neuronal_vs_glial.png`
- `output/figures/manuscript/fig_celltype_heatmap.png`
- `output/figures/manuscript/fig_celltype_by_direction.png`

**Key parameters:**
- Uses RNA assay (not SCT)
- Classification: >70% neuronal or >70% glial

**Dependencies:** 03, snRNA-seq data

---

### 05_snrnaseq/06b_score_snrnaseq.R

**Purpose:** Add module scores (Neuronal, Synaptic, Activity, Implant_Up, DAM, Complement) to the snRNA-seq dataset for downstream analysis.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- Module scores added in memory (no file output; informational script)

**Key parameters:**
- Uses RNA assay
- 6 signature modules defined
- Gene names converted to title case for rat

**Dependencies:** snRNA-seq data

---

### 05_snrnaseq/06b_stab_comparison.R

**Purpose:** Two-part analysis: (1) NPAS4 distance gradient from silicon timecourse, and (2) pseudobulk comparison of concordant genes in Implant vs Stab conditions using the snRNA-seq reference.

**Inputs:**
- `output/tables/timecourse/de_*.csv` (from 02b)
- `output/tables/comparison/validated_genes.csv` (from 03)
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/stab_comparison/npas4_distance_gradient.csv`
- `output/tables/stab_comparison/stab_vs_control_de.csv`
- `output/tables/stab_comparison/implant_vs_control_de.csv`
- `output/tables/stab_comparison/concordant_gene_classification.csv`
- `output/figures/manuscript/fig_npas4_distance_gradient.png`
- `output/figures/manuscript/fig_implant_vs_stab_scatter.png`
- `output/figures/manuscript/fig_gene_classification_bar.png`

**Key parameters:**
- DESeq2 pseudobulk by sample (orig.ident)
- Classification: Implant-Specific, Shared Injury, Stab-Only, Neither Significant

**Dependencies:** 02b, 03, snRNA-seq data

---

### 05_snrnaseq/06c_snrnaseq_stats.R

**Purpose:** Calculate summary statistics for the snRNA-seq dataset: cell counts by condition and cell type, percentage breakdowns.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/snrnaseq/snrnaseq_summary.csv`
- `output/tables/snrnaseq/snrnaseq_celltype_by_condition.csv`

**Dependencies:** snRNA-seq data

---

### 05_snrnaseq/06e_snrnaseq_pseudobulk.R

**Purpose:** Pseudobulk differential expression across all cell types using limma-voom. Runs Implant vs Control, Stab vs Control, and Implant vs Stab contrasts.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/snrnaseq/pseudobulk_DE_implant_vs_control.csv`
- `output/tables/snrnaseq/pseudobulk_DE_stab_vs_control.csv`
- `output/tables/snrnaseq/pseudobulk_DE_implant_vs_stab.csv`
- `output/tables/snrnaseq/pseudobulk_DE_summary.csv`

**Key parameters:**
- RNA assay, raw counts aggregated by orig.ident
- Significance: FDR < 0.05 AND |LFC| > 0.5
- limma-voom pipeline with edgeR normalization
- Minimum 10 cells per sample

**Dependencies:** snRNA-seq data

---

### 05_snrnaseq/06f_snrnaseq_signatures.R

**Purpose:** Comprehensive neuronal signature analysis. Scores 11 gene signature modules and tests for condition differences using Wilcoxon tests.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/snrnaseq/signature_scores_summary.csv`
- `output/tables/snrnaseq/signature_statistics.csv`
- `output/tables/snrnaseq/dysfunction_category_distribution.csv`
- `output/figures/snrnaseq/signature_heatmap.png`
- `output/figures/snrnaseq/signature_significant.png`
- `output/figures/snrnaseq/signature_violins.png`

**Key parameters:**
- 11 signatures: Activity_IEGs, Synaptic_Core, Catastrophic, Chronic, Preserved, Glutamatergic, GABAergic, Calcium_Signaling, Mitochondrial, Apoptosis, Survival
- Wilcoxon test with BH FDR correction
- Dysfunction categories: Catastrophic (<30%), Chronic (<70%), Preserved (>80%)

**Dependencies:** snRNA-seq data

---

### 05_snrnaseq/06g_microglial_states.R

**Purpose:** Microglial state analysis. Scores 4 state signatures (Homeostatic, Activated, DAM, Inflammatory) and assigns each microglia cell to its dominant state.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/snrnaseq/microglia_state_distribution.csv`
- `output/tables/snrnaseq/microglia_state_scores.csv`
- `output/tables/snrnaseq/microglia_state_statistics.csv`
- `output/figures/snrnaseq/microglia_state_distribution.png`
- `output/figures/snrnaseq/microglia_state_comparison.png`
- `output/figures/snrnaseq/microglia_state_violins.png`

**Key parameters:**
- State signatures: Homeostatic (P2ry12, Tmem119, ...), Activated (Cd68, Cd86, ...), DAM (Trem2, Apoe, ...), Inflammatory (Ccl2, Ccl5, ...)
- Dominant state: highest module score per cell (unless all scores < Q25)

**Dependencies:** snRNA-seq data

---

### 05_snrnaseq/06h_spp1_microglia_analysis.R

**Purpose:** Deep characterization of SPP1+ microglia. Identifies SPP1-high vs SPP1-low microglia, runs DE between them, performs GO/KEGG enrichment, and compares known microglial signatures.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/snrnaseq/spp1_microglia/spp1_by_condition.csv`
- `output/tables/snrnaseq/spp1_microglia/de_spp1high_vs_low.csv`
- `output/tables/snrnaseq/spp1_microglia/go_bp_spp1high_up.csv`
- `output/tables/snrnaseq/spp1_microglia/go_cc_spp1high_up.csv`
- `output/tables/snrnaseq/spp1_microglia/go_bp_spp1high_down.csv`
- `output/tables/snrnaseq/spp1_microglia/kegg_spp1high_up.csv`
- `output/tables/snrnaseq/spp1_microglia/signature_comparison.csv`
- `output/tables/snrnaseq/spp1_microglia/signature_statistics.csv`
- `output/tables/snrnaseq/spp1_microglia/spp1_high_by_condition.csv`
- `output/figures/snrnaseq/spp1_microglia/go_bp_spp1high_up.png`
- `output/figures/snrnaseq/spp1_microglia/signature_comparison_boxplot.png`

**Key parameters:**
- SPP1-high: Spp1 expression > 0 in RNA assay counts
- DE: FindAllMarkers (SPP1-high vs SPP1-low)
- GO enrichment: clusterProfiler enrichGO, org.Rn.eg.db

**Dependencies:** snRNA-seq data

---

### 05_snrnaseq/06i_complement_microglia_analysis.R

**Purpose:** Deep characterization of Complement+ (C1qa-high) microglia. Parallel analysis to 06h. Identifies complement-producing microglia and tests overlap with SPP1+ population.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/snrnaseq/complement_microglia/c1qa_by_condition.csv`
- `output/tables/snrnaseq/complement_microglia/de_comp_high_vs_low.csv`
- `output/tables/snrnaseq/complement_microglia/go_bp_comp_high_up.csv`
- `output/tables/snrnaseq/complement_microglia/go_bp_comp_high_down.csv`
- `output/tables/snrnaseq/complement_microglia/signature_comparison.csv`
- `output/tables/snrnaseq/complement_microglia/signature_statistics.csv`
- `output/tables/snrnaseq/complement_microglia/spp1_complement_overlap.csv`
- `output/figures/snrnaseq/complement_microglia/go_bp_comp_high_up.png`
- `output/figures/snrnaseq/complement_microglia/signature_comparison_boxplot.png`

**Key parameters:**
- Complement-high: top 25% C1qa expression (quartile-defined)
- Overlap test: Fisher's exact test for SPP1+ vs Complement+ independence

**Dependencies:** snRNA-seq data

---

### 05_snrnaseq/06j_microglia_populations_summary.R

**Purpose:** Comprehensive characterization of all microglial populations, including the 74% that are SPP1-/Complement-low. Scores 14 known microglial signatures.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/snrnaseq/microglia_populations/signature_by_population.csv`
- `output/tables/snrnaseq/microglia_populations/double_neg_dominant_by_condition.csv`
- `output/tables/snrnaseq/microglia_populations/population_by_condition.csv`
- `output/figures/snrnaseq/microglia_populations/signature_heatmap.png`
- `output/figures/snrnaseq/microglia_populations/population_by_condition.png`

**Key parameters:**
- 14 signatures: Homeostatic, DAM_Stage1, DAM_Stage2, Phagocytic, Inflammatory_M1, Anti_inflammatory_M2, Proliferative, IFN_Response, Lipid_Metabolism, Complement, Synapse_Pruning, SPP1_Module, Activated, Resting

**Dependencies:** snRNA-seq data

---

### 05_snrnaseq/06k_spatial_spp1_complement.R

**Purpose:** Spatial localization of SPP1+ vs Complement+ signatures in Visium data. Tests whether these two programs occupy different or overlapping spatial zones.

**Inputs:**
- `data/processed/spatial_seurat_list.RDS` (from 05a, or loads from raw)

**Outputs:**
- `output/tables/spatial/spp1_complement/spatial_correlation_summary.csv`
- `output/figures/manuscript/spp1_complement_spatial/chronic_no_stim_overview.png`
- `output/figures/manuscript/spp1_complement_spatial/zone_distribution_by_sample.png`
- `output/figures/manuscript/spp1_complement_spatial/spp1_complement_correlation_by_sample.png`

**Key parameters:**
- SPP1 genes: Spp1, Gpnmb, Lgals3, Fabp5, Igf1, Cd63, Lpl
- Complement genes: C1qa, C1qb, C1qc, C3
- Spatial correlation per sample

**Dependencies:** 05a

---

### 05_snrnaseq/06l_microglia_temporal_pseudobulk.R

**Purpose:** Temporal dynamics of SPP1+ and Complement+ microglial populations using pseudobulk (sample-level) statistics. Tests whether SPP1+ peaks early and declines.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/snrnaseq/microglia_temporal/` (sample summaries, temporal summaries)
- `output/figures/snrnaseq/microglia_temporal/spp1_temporal_pseudobulk.png`
- `output/figures/snrnaseq/microglia_temporal/complement_temporal_pseudobulk.png`
- `output/figures/snrnaseq/microglia_temporal/combined_temporal_pseudobulk.png`

**Key parameters:**
- Uses Duration/Timepoint metadata column
- Pseudobulk aggregation by sample (orig.ident)

**Dependencies:** snRNA-seq data

---

### 05_snrnaseq/06l_spp1_implant_vs_stab.R

**Purpose:** Direct pseudobulk DE test of SPP1 expression in microglia comparing Implant vs Stab conditions. Formal test for implant-specificity.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/snrnaseq/spp1_specificity/spp1_pairwise_comparisons.csv`
- `output/tables/snrnaseq/spp1_specificity/implant_vs_stab_de.csv`

**Key parameters:**
- Microglia only (celltype_l1 == "Microglia")
- RNA assay, pseudobulk by sample
- limma-voom with edgeR normalization

**Dependencies:** snRNA-seq data

---

### 05_snrnaseq/06m_spp1_proportion_test.R

**Purpose:** Test SPP1+ cell proportions across conditions (Control, Implant, Stab) at the sample level. For rare markers, proportion of positive cells is the appropriate metric.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/snrnaseq/spp1_specificity/spp1_proportion_tests.csv`
- `output/tables/snrnaseq/spp1_specificity/spp1_by_sample.csv`

**Key parameters:**
- SPP1+: Spp1 count > 0 in RNA assay
- Sample-level proportion test

**Dependencies:** snRNA-seq data

---

### 05_snrnaseq/06n_neuronal_pseudobulk_de.R

**Purpose:** Pseudobulk differential expression for neurons only. Identifies neuron-specific gene changes at the electrode interface.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/snrnaseq/neurons/de_implant_vs_control.csv`
- `output/tables/snrnaseq/neurons/key_neuronal_genes.csv`
- `output/tables/snrnaseq/neurons/go_bp_down.csv`
- `output/tables/snrnaseq/neurons/go_bp_up.csv`

**Key parameters:**
- Neuron subset: celltype_l1 in (Neuron_Excitatory, Neuron_Inhibitory, etc.)
- limma-voom pseudobulk
- GO enrichment of up/down regulated genes

**Dependencies:** snRNA-seq data

---

### 05_snrnaseq/06n_neuronal_subtype_de.R

**Purpose:** Pseudobulk DE for each neuronal subtype separately (Excitatory_Projection, Inhibitory_PV, Inhibitory_CCK). Tests differential vulnerability.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/snrnaseq/neurons/de_{subtype}_implant_vs_ctrl.csv` (3 files)
- `output/tables/snrnaseq/neurons/key_genes_by_subtype.csv`

**Key parameters:**
- Uses celltype_l3 annotations
- Subtypes: Neuron_Excitatory_Projection (44,807 cells), Neuron_Inhibitory_PV (7,704), Neuron_Inhibitory_CCK (6,634)
- Cell counts: ~45K excitatory, ~7K PV+, ~7K CCK+ (asymmetric representation)
- Minimum 30 cells per sample for pseudobulk

**Dependencies:** snRNA-seq data

---

### 05_snrnaseq/06o_oligo_subtype_de.R

**Purpose:** Pseudobulk DE for each oligodendrocyte subtype separately (Myelinating, Precursor/OPC, Premyelinating).

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/snrnaseq/oligodendrocytes/de_{subtype}_implant_vs_ctrl.csv` (3 files)
- `output/tables/snrnaseq/oligodendrocytes/key_genes_by_subtype.csv`

**Key parameters:**
- Uses celltype_l2 annotations
- Subtypes: Myelinating (8,791 cells), Precursor (4,972), Premyelinating (659)
- Minimum 10 cells per sample for pseudobulk

**Dependencies:** snRNA-seq data

---

### 05_snrnaseq/06o_oligodendrocyte_pseudobulk_de.R

**Purpose:** Pseudobulk DE for all oligodendrocyte lineage cells combined, with GO enrichment of dysregulated genes.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/snrnaseq/oligodendrocytes/de_implant_vs_control.csv`
- GO enrichment results in `output/tables/snrnaseq/oligodendrocytes/`

**Key parameters:**
- All oligodendrocyte lineage cells (Myelinating + Precursor + Premyelinating)
- limma-voom pseudobulk
- GO enrichment via clusterProfiler

**Dependencies:** snRNA-seq data

---

## Stage 08: Astrocyte Analysis (`R/06_astrocyte/`)

---

### 06_astrocyte/08_astrocyte_analysis.R

**Purpose:** Comprehensive astrocyte analysis addressing asymmetric focus on microglia. Covers reactive vs homeostatic states, complement (C3) production by astrocytes, glial scar markers, and astrocyte-specific pseudobulk DE.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/snrnaseq/astrocytes/state_scores_by_condition.csv`
- `output/tables/snrnaseq/astrocytes/state_statistics.csv`
- `output/tables/snrnaseq/astrocytes/marker_expression.csv`
- `output/tables/snrnaseq/astrocytes/de_implant_vs_control.csv`
- `output/tables/snrnaseq/astrocytes/de_stab_vs_control.csv`
- `output/tables/snrnaseq/astrocytes/c3_astrocyte_vs_microglia.csv`
- `output/tables/snrnaseq/astrocytes/c1q_astrocyte_vs_microglia.csv`
- `output/tables/snrnaseq/astrocytes/summary_stats.rds`
- `output/figures/snrnaseq/astrocytes/state_violins.png`
- `output/figures/snrnaseq/astrocytes/c3_comparison.png`

**Key parameters:**
- Astrocyte subset: celltype_l1 == "Astrocyte"
- State signatures: A1 (neurotoxic), A2 (neuroprotective), Pan-reactive, Homeostatic
- Key comparison: C3 production in astrocytes vs microglia

**Dependencies:** snRNA-seq data

---

### 06_astrocyte/08b_astrocyte_spatial.R

**Purpose:** Spatial analysis of astrocyte signatures and C3 expression. Tests the two-cell complement model where microglia produce C1q and astrocytes produce C3.

**Inputs:**
- Visium data from `data/external/stRNAseq/` (loaded from source)

**Outputs:**
- `output/tables/spatial/astrocyte/spatial_correlation_summary.csv`
- `output/tables/spatial/astrocyte/two_cell_model_summary.csv`
- `output/figures/manuscript/astrocyte_spatial/two_cell_model_spatial.png`
- `output/figures/manuscript/astrocyte_spatial/two_cell_model_correlations.png`
- `output/figures/manuscript/astrocyte_spatial/two_cell_model_forest.png`

**Key parameters:**
- Microglia C1q genes: C1qa, C1qb, C1qc
- Astrocyte C3: C3
- 6 Visium sections

**Dependencies:** Visium raw data

---

## Stage 09: Robustness Analyses (`R/07_robustness/`)

---

### 07_robustness/09a_spatial_density_permutation.R

**Purpose:** Density-controlled permutation test for spatial colocalization of SPP1 and Complement signatures. Controls for UMI density gradients in the null model.

**Inputs:**
- `data/processed/spatial_scored.RDS` (from 05b) or `spatial_seurat_list.RDS` (from 05a)

**Outputs:**
- `output/tables/spatial/permutation/density_permutation_results.csv`
- `output/tables/spatial/permutation/null_distributions.csv`
- `output/figures/spatial/permutation/` (null distribution plots)

**Key parameters:**
- SPP1 genes: Spp1, Gpnmb, Lgals3, Fabp5, Igf1, Cd63, Lpl
- Complement genes: C1qa, C1qb, C1qc, C3
- Null model: permute signature scores within UMI-count bins (preserves density structure)

**Dependencies:** 05a or 05b

---

### 07_robustness/09b_killzone_dimensions_robust.R

**Purpose:** Robust kill zone dimension estimation with Otsu thresholding and bootstrap confidence intervals. Uses data-driven thresholding as an alternative to quartile-based boundaries.

**Inputs:**
- `data/processed/spatial_scored.RDS` (from 05b) or `spatial_seurat_list.RDS` (from 05a)

**Outputs:**
- `output/tables/spatial/dimensions/killzone_dimensions_robust.csv`
- `output/tables/spatial/dimensions/condition_summary.csv`
- `output/tables/spatial/dimensions/expansion_test.csv`
- `output/tables/spatial/dimensions/bootstrap_distributions.csv`
- `output/figures/spatial/dimensions/` (radii plots, expansion test)
- `output/figures/manuscript/killzone/` (publication figures)

**Key parameters:**
- Otsu threshold: automatic, data-driven kill zone boundary
- Bootstrap: resampling for CI estimation

**Dependencies:** 05a or 05b

---

### 07_robustness/09c_population_bootstrap.R

**Purpose:** Bootstrap confidence intervals for microglial population overlap statistics, particularly the SPP1+ vs Complement+ odds ratio and population proportions.

**Inputs:**
- snRNA-seq Seurat object (via SNRNASEQ_PATH)

**Outputs:**
- `output/tables/snrnaseq/population_bootstrap/summary.csv`
- `output/tables/snrnaseq/population_bootstrap/boot_distributions.csv`
- `output/tables/snrnaseq/population_bootstrap/stability.csv`
- `output/tables/snrnaseq/population_bootstrap/condition_results.csv`
- `output/figures/snrnaseq/population_bootstrap/` (bootstrap distribution plots)

**Key parameters:**
- SPP1+: Spp1 > 0 in RNA counts
- Complement-high: top 25% C1qa expression
- Bootstrap iterations: 1000

**Dependencies:** snRNA-seq data

---

## Stage 10: Sensitivity Analyses (`R/07_robustness/`)

---

### 07_robustness/10a_concordance_sensitivity.R

**Purpose:** Sensitivity analysis testing how the number of concordant genes changes at relaxed significance thresholds. Determines whether 50 genes is a floor or ceiling.

**Inputs:**
- `output/tables/comparison/cross_platform_validation.csv` (from 03a)

**Outputs:**
- `output/tables/comparison/concordance_sensitivity.csv`
- `output/tables/comparison/concordance_relaxed_genes.csv` (if relaxed set generated)

**Key parameters:**
- Tiers: FDR<0.05 both (strict), FDR<0.05/FDR<0.10, FDR<0.05/nominal p<0.05, direction-only
- 557 matched genes as base set

**Dependencies:** 03a

---

### 07_robustness/10b_bootstrap_ci_check.R

**Purpose:** Verification of bootstrap CI methods for Table 1 neuronal suppression values. Compares three CI estimation methods for robustness.

**Inputs:**
- `output/tables/snrnaseq/signature_scores_summary.csv` (from 06f)

**Outputs:**
- `output/tables/snrnaseq/bootstrap_ci_comparison.csv`

**Key parameters:**
- 5 key programs: Calcium_Signaling, Survival, Synaptic_Core, Glutamatergic, Activity_IEGs
- Methods compared: Delta (normal approximation), Percentile bootstrap, BCa bootstrap
- Compares three CI estimation methods for robustness

**Dependencies:** 06f
