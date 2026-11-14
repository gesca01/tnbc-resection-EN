## Collects the per-mouse scores and summary numbers for Fig S3 into one RDS.
## Uses the 15-gene day-7 treated signature with the treatment-arm offset removed
## (see _day7_signature_common.R).

options(stringsAsFactors = FALSE)
.libPaths(c("./NCT_paper_pipeline/renv_lib", .libPaths()))
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(edgeR); library(limma)
  library(matrixStats); library(AnnotationDbi); library(org.Mm.eg.db)
  library(pROC); library(cluster)
})
set.seed(123)

PROJ <- "."
SRC  <- file.path(PROJ, "Day7_Treated_Signature_Applied_to_Untreated")
OUT  <- file.path(PROJ, "Supplementary_Figure_Day7_Transfer")
dir.create(file.path(OUT, "data"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "tables"), recursive = TRUE, showWarnings = FALSE)

source(file.path(SRC, "scripts", "_day7_signature_common.R"))

## ---- signature space (all 23 Day 7 mice, arm-centered) ----
S    <- build_day7_signature_space(PROJ)
M    <- S$M                       # 15 genes x 23 mice, arm-centered z-scores
meta <- S$meta
meta$eta <- S$eta                 # risk score: higher = more like a treated non-survivor
meta$Arm    <- factor(meta$TxLabel,  levels = c("Treated", "Untreated"))
meta$Actual <- factor(meta$GroupBin, levels = c("Survivor", "NonSurvivor"))
stopifnot(nrow(meta) == 23, length(S$sig) == 15)
message("Signature genes: ", paste(S$sig, collapse = ", "))

is_ut <- meta$Arm == "Untreated"

## ---- risk-axis call: median split WITHIN the untreated cohort ----
## Untreated mice cannot sit on an absolute treated scale - they were never
## treated - so the call is a ranking among themselves, cut at their median.
med_ut <- median(meta$eta[is_ut])
meta$Call <- NA_character_
meta$Call[is_ut] <- ifelse(meta$eta[is_ut] >= med_ut,
                           "Predicted non-survivor", "Predicted survivor")
meta$Call <- factor(meta$Call, levels = c("Predicted survivor", "Predicted non-survivor"))
message("Median cut = ", round(med_ut, 3), "  -> ",
        sum(meta$Call == "Predicted survivor", na.rm = TRUE), " pred-survivor / ",
        sum(meta$Call == "Predicted non-survivor", na.rm = TRUE), " pred-non-survivor")

## ---- performance of the score against actual outcome ----
auc_of <- function(y, s) {
  r  <- suppressWarnings(roc(y, s, levels = c(0, 1), direction = "<", quiet = TRUE))
  ci <- suppressWarnings(ci.auc(r, method = "delong"))
  c(auc = as.numeric(auc(r)), lo = as.numeric(ci)[1], hi = as.numeric(ci)[3])
}
y_ut <- as.integer(meta$Actual[is_ut]  == "NonSurvivor")
y_tr <- as.integer(meta$Actual[!is_ut] == "NonSurvivor")
auc_ut <- auc_of(y_ut, meta$eta[is_ut])
auc_tr <- auc_of(y_tr, meta$eta[!is_ut])
w_ut   <- suppressWarnings(wilcox.test(meta$eta[is_ut] ~ y_ut))$p.value
concord <- mean((meta$Call[is_ut] == "Predicted non-survivor") == (y_ut == 1))

perf <- data.frame(
  cohort  = c("Treated (training, refit-on-all)", "Untreated (out-of-cohort)"),
  n       = c(sum(!is_ut), sum(is_ut)),
  AUC     = c(auc_tr["auc"], auc_ut["auc"]),
  CI_low  = c(auc_tr["lo"],  auc_ut["lo"]),
  CI_high = c(auc_tr["hi"],  auc_ut["hi"]),
  wilcox_p = c(NA, w_ut))
print(perf, digits = 3)
message("Call vs actual untreated outcome: ", round(concord, 3),
        "  (", sum((meta$Call[is_ut] == "Predicted non-survivor") == (y_ut == 1)), "/11)")

## ---- unsupervised structure: label-blind Ward.D2 + silhouette ----
## Nothing about outcome, arm or the score is consulted until after clusters exist.
cluster_block <- function(Msub, ks = 2:6) {
  dm <- dist(t(Msub)); hc <- hclust(dm, method = "ward.D2")
  sil <- vapply(ks, function(k) mean(silhouette(cutree(hc, k), dm)[, "sil_width"]), numeric(1))
  best <- ks[which.max(sil)]
  list(hc = hc, order = hc$order, sil_tbl = data.frame(k = ks, sil = sil),
       best_k = best, cl = factor(paste0("C", cutree(hc, best))))
}
cl23 <- cluster_block(M)
cl_ut <- cluster_block(M[, is_ut, drop = FALSE])
meta$Cluster23 <- cl23$cl
message("all 23: k=", cl23$best_k, " sil=", round(max(cl23$sil_tbl$sil), 3),
        " | untreated: k=", cl_ut$best_k, " sil=", round(max(cl_ut$sil_tbl$sil), 3))

## ---- PCA of the same 15-gene space; PC1 oriented to run with the risk score ----
pc  <- prcomp(t(M), scale. = FALSE)
varx <- 100 * pc$sdev^2 / sum(pc$sdev^2)
if (cor(pc$x[, 1], meta$eta) < 0) pc$x[, 1] <- -pc$x[, 1]   # PC1 sign is arbitrary
meta$PC1 <- pc$x[, 1]; meta$PC2 <- pc$x[, 2]
message(sprintf("PC1 %.1f%% of variance; cor(PC1, risk score) = %.3f",
                varx[1], cor(meta$PC1, meta$eta)))

## ---- gene coefficients, for the heatmap row annotation ----
genes <- data.frame(gene = S$sig, coef = S$beta,
                    direction = ifelse(S$beta > 0, "NS-high", "S-high"))

## ---- write ----
write_csv(meta %>%
            dplyr::select(Core_Sample_Id, Sample_Name, Arm, Actual, eta, PC1, PC2, Call, Cluster23) %>%
            arrange(Arm, desc(eta)),
          file.path(OUT, "tables", "S_per_mouse_scores_and_calls.csv"))
write_csv(perf,  file.path(OUT, "tables", "S_transfer_performance.csv"))
write_csv(genes, file.path(OUT, "tables", "S_signature_genes.csv"))
write_csv(bind_rows(cl23$sil_tbl %>% mutate(cohort = "All 23"),
                    cl_ut$sil_tbl %>% mutate(cohort = "Untreated (n=11)")),
          file.path(OUT, "tables", "S_silhouette_by_k.csv"))

saveRDS(list(M = M, meta = meta, genes = genes, perf = perf, med_ut = med_ut,
             concord = concord, wilcox_p_ut = w_ut,
             cl23 = cl23, cl_ut = cl_ut, varx = varx),
        file.path(OUT, "data", "panel_data.rds"))
message("Wrote data/panel_data.rds")
