# 06j_microglia_populations_summary.R
# Comprehensive characterization of all microglial populations
# What are the 74% that are SPP1-/Complement-Low?

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/config.R")

cat("=== Complete Microglial Population Characterization ===\n\n")

dir.create(file.path(OUT_TABLES_SNRNASEQ, "microglia_populations"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_FIGURES_SNRNASEQ, "microglia_populations"), recursive = TRUE, showWarnings = FALSE)

# Load snRNA-seq data
SNRNASEQ_SOURCE <- SNRNASEQ_PATH  # defined in config.R
cat("Loading snRNA-seq data...\n")
seu <- readRDS(SNRNASEQ_SOURCE)
DefaultAssay(seu) <- "RNA"

# Subset to microglia
celltype_col <- if ("celltype_l2" %in% colnames(seu@meta.data)) "celltype_l2" else "celltype_l1"
microglia_types <- unique(seu@meta.data[[celltype_col]])[grepl("Microglia|microglia|MG", unique(seu@meta.data[[celltype_col]]), ignore.case = TRUE)]
seu_mg <- subset(seu, cells = colnames(seu)[seu@meta.data[[celltype_col]] %in% microglia_types])
# CRITICAL: Explicitly set RNA assay after subset for accurate gene detection
DefaultAssay(seu_mg) <- "RNA"
cat(sprintf("Microglia: %d cells (using %s assay)\n", ncol(seu_mg), DefaultAssay(seu_mg)))

# Define all signatures
microglial_signatures <- list(
  Homeostatic = c("P2ry12", "Tmem119", "Cx3cr1", "Siglech", "Hexb", "Cst3", "Sall1"),
  DAM_Stage1 = c("Tyrobp", "Ctsd", "Lyz2", "Ctsb", "Apoe"),
  DAM_Stage2 = c("Trem2", "Axl", "Cst7", "Lpl", "Cd9", "Csf1", "Itgax", "Clec7a", "Lilrb4"),
  Phagocytic = c("Cd68", "Mertk", "Axl", "Tyro3", "Megf10", "Gulp1"),
  Inflammatory_M1 = c("Il1b", "Tnf", "Il6", "Nos2", "Cd86", "Ccl2", "Cxcl10"),
  Anti_inflammatory_M2 = c("Arg1", "Mrc1", "Cd163", "Il10", "Tgfb1", "Chil3"),
  Proliferative = c("Mki67", "Top2a", "Pcna", "Cdk1", "Ccnb1"),
  IFN_Response = c("Ifit1", "Ifit2", "Ifit3", "Isg15", "Mx1", "Stat1", "Irf7"),
  Lipid_Metabolism = c("Apoe", "Lpl", "Fabp5", "Abca1", "Abcg1", "Plin2"),
  Complement = c("C1qa", "C1qb", "C1qc", "C3"),
  Synapse_Pruning = c("C1qa", "C1qb", "C1qc", "C3", "Trem2", "Tyrobp"),
  SPP1_Module = c("Spp1", "Gpnmb", "Lgals3", "Fabp5", "Igf1", "Cd63", "Lpl"),
  Activated = c("Cd68", "Cd86", "Il1b", "Tnf", "Nos2", "Ptgs2"),
  Resting = c("P2ry12", "Tmem119", "Cx3cr1", "Hexb", "Cst3")
)

# Score all signatures
available_genes <- rownames(seu_mg)
for (sig_name in names(microglial_signatures)) {
  genes <- microglial_signatures[[sig_name]]
  genes_title <- tools::toTitleCase(tolower(genes))
  present <- intersect(genes_title, available_genes)

  if (length(present) >= 3) {
    seu_mg <- AddModuleScore(seu_mg, features = list(present), name = sig_name, seed = 42)
    colnames(seu_mg@meta.data)[ncol(seu_mg@meta.data)] <- sig_name
  }
}

# Get SPP1 and complement expression
# NOTE: Use RNA assay explicitly - SCT may distort sparse gene detection
spp1_expr <- GetAssayData(seu_mg, assay = "RNA", layer = "data")["Spp1", ]
seu_mg$spp1_pos <- spp1_expr > 0

# Define complement groups
comp_q75 <- quantile(seu_mg$Complement, 0.75)
comp_q25 <- quantile(seu_mg$Complement, 0.25)
seu_mg$comp_group <- case_when(
  seu_mg$Complement >= comp_q75 ~ "Comp_High",
  seu_mg$Complement <= comp_q25 ~ "Comp_Low",
  TRUE ~ "Comp_Mid"
)

# Create 4-population classification
seu_mg$population <- case_when(
  seu_mg$spp1_pos & seu_mg$comp_group == "Comp_High" ~ "SPP1+/Comp+",
  seu_mg$spp1_pos & seu_mg$comp_group != "Comp_High" ~ "SPP1+/Comp-",
  !seu_mg$spp1_pos & seu_mg$comp_group == "Comp_High" ~ "SPP1-/Comp+",
  TRUE ~ "SPP1-/Comp-"
)

cat("\n=== POPULATION DISTRIBUTION ===\n")
pop_dist <- table(seu_mg$population)
print(pop_dist)
cat(sprintf("\nPercentages:\n"))
print(round(100 * prop.table(pop_dist), 1))

# ============================================================
# Characterize each population by signatures
# ============================================================

cat("\n=== SIGNATURE SCORES BY POPULATION ===\n")

sig_cols <- intersect(names(microglial_signatures), colnames(seu_mg@meta.data))

pop_signatures <- seu_mg@meta.data %>%
  group_by(population) %>%
  summarize(
    n_cells = n(),
    across(all_of(sig_cols), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )

cat("\nMean signature scores by population:\n")
print(pop_signatures)

# Transpose for easier reading
pop_sig_long <- pop_signatures %>%
  pivot_longer(-c(population, n_cells), names_to = "signature", values_to = "score") %>%
  pivot_wider(names_from = population, values_from = score)

write.csv(pop_sig_long, file.path(OUT_TABLES_SNRNASEQ, "microglia_populations", "signature_by_population.csv"), row.names = FALSE)

cat("\nSignature scores (wide format):\n")
print(pop_sig_long)

# ============================================================
# What dominates the SPP1-/Comp- population?
# ============================================================

cat("\n=== SPP1-/COMP- POPULATION ANALYSIS (74% of microglia) ===\n")

double_neg <- seu_mg@meta.data %>% filter(population == "SPP1-/Comp-")
cat(sprintf("SPP1-/Comp- cells: %d (%.1f%%)\n", nrow(double_neg), 100 * nrow(double_neg) / ncol(seu_mg)))

# Find dominant signature for each cell
sig_scores <- double_neg[, sig_cols]
double_neg$dominant_signature <- apply(sig_scores, 1, function(x) names(x)[which.max(x)])

cat("\nDominant signature in SPP1-/Comp- population:\n")
dom_dist <- table(double_neg$dominant_signature)
print(sort(dom_dist, decreasing = TRUE))
cat("\nPercentages:\n")
print(round(100 * sort(prop.table(dom_dist), decreasing = TRUE), 1))

# By condition
cat("\nSPP1-/Comp- by condition:\n")
print(table(double_neg$Condition))

# Dominant signature by condition in SPP1-/Comp-
dom_by_cond <- double_neg %>%
  group_by(Condition, dominant_signature) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  arrange(Condition, desc(pct))

cat("\nDominant signature by condition (SPP1-/Comp-):\n")
print(dom_by_cond)

write.csv(dom_by_cond, file.path(OUT_TABLES_SNRNASEQ, "microglia_populations", "double_neg_dominant_by_condition.csv"), row.names = FALSE)

# ============================================================
# Compare all 4 populations
# ============================================================

cat("\n=== ALL POPULATION COMPARISON ===\n")

# Dominant signature for ALL cells
all_sig_scores <- seu_mg@meta.data[, sig_cols]
seu_mg$dominant_signature <- apply(all_sig_scores, 1, function(x) names(x)[which.max(x)])

# Cross-tabulation
cat("\nDominant signature by population:\n")
cross_tab <- table(seu_mg$dominant_signature, seu_mg$population)
print(cross_tab)

cat("\nPercentages within each population:\n")
print(round(100 * prop.table(cross_tab, margin = 2), 1))

# By condition
cat("\n\nPopulation distribution by condition:\n")
pop_by_cond <- seu_mg@meta.data %>%
  group_by(Condition, population) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(pct = 100 * n / sum(n))

print(pop_by_cond %>% arrange(Condition, desc(pct)))

write.csv(pop_by_cond, file.path(OUT_TABLES_SNRNASEQ, "microglia_populations", "population_by_condition.csv"), row.names = FALSE)

# ============================================================
# Visualization
# ============================================================

# Heatmap of signatures by population
heatmap_data <- pop_sig_long %>%
  pivot_longer(-signature, names_to = "population", values_to = "score") %>%
  filter(population != "n_cells")

p_heat <- ggplot(heatmap_data, aes(x = population, y = signature, fill = score)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(score, 2)), size = 3) +
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C", midpoint = 0) +
  labs(title = "Microglial Signatures by Population",
       x = NULL, y = NULL, fill = "Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_FIGURES_SNRNASEQ, "microglia_populations", "signature_heatmap.png"),
       p_heat, width = 10, height = 10)

# Bar plot of population by condition
p_pop <- ggplot(pop_by_cond, aes(x = Condition, y = pct, fill = population)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "SPP1-/Comp-" = "#BDBDBD",
    "SPP1-/Comp+" = "#984EA3",
    "SPP1+/Comp-" = "#E41A1C",
    "SPP1+/Comp+" = "#FF7F00"
  )) +
  labs(title = "Microglial Population Distribution by Condition",
       x = NULL, y = "% of Microglia", fill = "Population") +
  theme_minimal()

ggsave(file.path(OUT_FIGURES_SNRNASEQ, "microglia_populations", "population_by_condition.png"),
       p_pop, width = 8, height = 6)

# ============================================================
# SUMMARY
# ============================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("COMPLETE MICROGLIAL POPULATION SUMMARY\n")
cat(strrep("=", 70), "\n")

cat("\n[1] POPULATION DISTRIBUTION:\n")
for (pop in names(sort(pop_dist, decreasing = TRUE))) {
  cat(sprintf("    %s: %d cells (%.1f%%)\n", pop, pop_dist[pop], 100 * pop_dist[pop] / sum(pop_dist)))
}

cat("\n[2] SPP1-/COMP- POPULATION (The Majority):\n")
cat(sprintf("    Total: %d cells (%.1f%% of all microglia)\n", nrow(double_neg), 100 * nrow(double_neg) / ncol(seu_mg)))
cat("    Dominant signatures:\n")
top_dom <- head(sort(dom_dist, decreasing = TRUE), 5)
for (i in seq_along(top_dom)) {
  cat(sprintf("      %d. %s: %.1f%%\n", i, names(top_dom)[i], 100 * top_dom[i] / sum(dom_dist)))
}

cat("\n[3] KEY INSIGHT:\n")
cat("    The 'rest' of microglia (SPP1-/Comp-) are predominantly:\n")
homeo_pct <- 100 * sum(double_neg$dominant_signature == "Homeostatic") / nrow(double_neg)
resting_pct <- 100 * sum(double_neg$dominant_signature == "Resting") / nrow(double_neg)
cat(sprintf("    - Homeostatic: %.1f%%\n", homeo_pct))
cat(sprintf("    - Resting: %.1f%%\n", resting_pct))
cat(sprintf("    - Combined quiescent: %.1f%%\n", homeo_pct + resting_pct))

cat("\n[4] BIOLOGICAL INTERPRETATION:\n")
cat("    The majority of microglia remain in a quiescent/homeostatic state.\n")
cat("    Only ~25% become Complement+ (synapse pruning) and ~1% become SPP1+.\n")
cat("    This suggests the FBR involves specialized subpopulations, not\n")
cat("    wholesale microglial activation.\n")

cat("\nOutputs saved to:\n")
cat(sprintf("  %s/microglia_populations/\n", OUT_TABLES_SNRNASEQ))
cat(sprintf("  %s/microglia_populations/\n", OUT_FIGURES_SNRNASEQ))
