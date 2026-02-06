# 02_polyimide_limma.R
# limma analysis: Implant vs Control (primary) + three-way comparison
# Output: output/tables/polyimide_limma_results.csv, implant_specific.csv, shared_injury.csv

library(oligo)
library(limma)
library(dplyr)
library(affycoretools)
library(pd.clariom.s.rat)
source("R/config.R")

# Load and normalize
cel_files <- list.files("data/Internal/Polyimide_Microarray", pattern = "\\.CEL$", full.names = TRUE)
raw <- read.celfiles(cel_files)
eset <- rma(raw)

# Annotate with gene symbols
eset <- annotateEset(eset, pd.clariom.s.rat)
fdata <- fData(eset)
cat(sprintf("Annotated: %d probes, %d with symbols\n", nrow(fdata), sum(!is.na(fdata$SYMBOL))))

# Aggregate to gene level (average probes per gene)
expr <- exprs(eset)
symbols <- fdata$SYMBOL
valid <- !is.na(symbols) & symbols != ""
expr <- expr[valid, ]
symbols <- symbols[valid]

# Average by gene
expr_df <- data.frame(SYMBOL = symbols, expr, check.names = FALSE)
expr_agg <- aggregate(. ~ SYMBOL, data = expr_df, FUN = mean)
rownames(expr_agg) <- expr_agg$SYMBOL
expr_mat <- as.matrix(expr_agg[, -1])
cat(sprintf("Genes after aggregation: %d\n", nrow(expr_mat)))

# Sample metadata
samples <- colnames(expr_mat)
condition <- case_when(
  grepl("Control", samples) ~ "Control",
  grepl("Week1", samples) ~ "Implant",
  grepl("Stab", samples) ~ "Stab"
)

# Primary: Implant vs Control
idx <- condition %in% c("Control", "Implant")
design <- model.matrix(~ factor(condition[idx], levels = c("Control", "Implant")))
fit <- lmFit(expr_mat[, idx], design)
fit <- eBayes(fit)
res <- topTable(fit, coef = 2, number = Inf)
res$gene <- rownames(res)
res$significant <- res$adj.P.Val < FDR_THRESH
res$direction <- ifelse(res$logFC > 0 & res$significant, "Up",
                        ifelse(res$logFC < 0 & res$significant, "Down", "NS"))

cat(sprintf("Implant vs Control DEGs: %d (Up: %d, Down: %d)\n",
            sum(res$significant), sum(res$direction == "Up"), sum(res$direction == "Down")))
write.csv(res, file.path(OUT_TABLES_DEG, "polyimide_limma_results.csv"), row.names = FALSE)

# Three-way comparison
cond_f <- factor(condition, levels = c("Control", "Implant", "Stab"))
design_all <- model.matrix(~ 0 + cond_f)
colnames(design_all) <- c("Control", "Implant", "Stab")
fit_all <- lmFit(expr_mat, design_all)
contrast_mat <- makeContrasts(
  Implant_vs_Control = Implant - Control,
  Stab_vs_Control = Stab - Control,
  Implant_vs_Stab = Implant - Stab,
  levels = design_all
)
fit_con <- contrasts.fit(fit_all, contrast_mat)
fit_con <- eBayes(fit_con)

res_ic <- topTable(fit_con, coef = "Implant_vs_Control", number = Inf)
res_sc <- topTable(fit_con, coef = "Stab_vs_Control", number = Inf)
res_is <- topTable(fit_con, coef = "Implant_vs_Stab", number = Inf)

genes <- rownames(res_ic)
three_way <- data.frame(
  gene = genes,
  lfc_impl_ctrl = res_ic[genes, "logFC"],
  fdr_impl_ctrl = res_ic[genes, "adj.P.Val"],
  lfc_stab_ctrl = res_sc[genes, "logFC"],
  fdr_stab_ctrl = res_sc[genes, "adj.P.Val"],
  lfc_impl_stab = res_is[genes, "logFC"],
  fdr_impl_stab = res_is[genes, "adj.P.Val"]
)

three_way$sig_ic <- three_way$fdr_impl_ctrl < FDR_THRESH
three_way$sig_sc <- three_way$fdr_stab_ctrl < FDR_THRESH
three_way$sig_is <- three_way$fdr_impl_stab < FDR_THRESH

three_way$category <- with(three_way, case_when(
  sig_ic & (sig_is | !sig_sc) & lfc_impl_ctrl > 0 ~ "Implant_Up",
  sig_ic & (sig_is | !sig_sc) & lfc_impl_ctrl < 0 ~ "Implant_Down",
  sig_ic & sig_sc & sign(lfc_impl_ctrl) == sign(lfc_stab_ctrl) & !sig_is & lfc_impl_ctrl > 0 ~ "Shared_Up",
  sig_ic & sig_sc & sign(lfc_impl_ctrl) == sign(lfc_stab_ctrl) & !sig_is & lfc_impl_ctrl < 0 ~ "Shared_Down",
  TRUE ~ "NS"
))

implant_specific <- three_way[three_way$category %in% c("Implant_Up", "Implant_Down"), ]
shared_injury <- three_way[three_way$category %in% c("Shared_Up", "Shared_Down"), ]

cat(sprintf("\nImplant-specific: %d (Up: %d, Down: %d)\n",
            nrow(implant_specific), sum(implant_specific$category == "Implant_Up"),
            sum(implant_specific$category == "Implant_Down")))
cat(sprintf("Shared injury: %d (Up: %d, Down: %d)\n",
            nrow(shared_injury), sum(shared_injury$category == "Shared_Up"),
            sum(shared_injury$category == "Shared_Down")))

write.csv(implant_specific, file.path(OUT_TABLES_COMPARISON, "implant_specific.csv"), row.names = FALSE)
write.csv(shared_injury, file.path(OUT_TABLES_COMPARISON, "shared_injury.csv"), row.names = FALSE)
write.csv(three_way, file.path(OUT_TABLES_COMPARISON, "polyimide_three_way.csv"), row.names = FALSE)

cat("\nSaved: output/tables/deg/ and output/tables/comparison/\n")
