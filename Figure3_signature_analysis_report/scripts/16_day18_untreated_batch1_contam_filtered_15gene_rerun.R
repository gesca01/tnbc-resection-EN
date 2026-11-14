## Day-18 untreated signature (Fig S2, Table S7), sequencing batch 1 only.
## Blood/erythroid and epithelial/keratinocyte genes (likely tissue contamination)
## are removed before ranking, and the signature size is limited to 12-18 genes.
## Same pipeline as the main signature run otherwise.

options(stringsAsFactors = FALSE)
set.seed(123)

PROJ <- "."
PIPE <- file.path(PROJ, "NCT_paper_pipeline")

.libPaths(c(file.path(PIPE, "renv_lib"), .libPaths()))

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(matrixStats); library(limma)
  library(glmnet); library(pROC); library(edgeR)
})

INPUT_COUNTS <- file.path(PIPE, "01_inputs", "merged_counts.tsv")
INPUT_META   <- file.path(PIPE, "01_inputs", "merged_metadata_all.tsv")
SIG_ROOT       <- file.path(PIPE, "02_cache", "elastic_net_day7_day18_source_match",
                            "results_EN_signatures")
BATCH1_ROOT    <- file.path(PIPE, "02_cache", "elastic_net_day7_day18_source_match",
                            "results_EN_day18_untreated_batch1")
OUTDIR         <- file.path(PIPE, "02_cache", "elastic_net_day7_day18_source_match",
                            "results_EN_day18_untreated_batch1_contam_filtered_15gene",
                            "Day18_Untreated_SvNS")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

K_GRID  <- c(5, 7, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200)
ALPHAS  <- c(0.01, 0.02, 0.03, seq(0.05, 1.0, by = 0.05))
SIG_MIN <- 12L   # narrowed from the default (8-30) to target ~15 genes
SIG_MAX <- 18L
RKFOLD_K       <- 5L
RKFOLD_REPEATS <- 10L
N_PERM <- 500L

TARGET_BATCH <- "batch_1"

## ---- contamination gene families to exclude before ranking ----
## Same blood/erythroid family as script 13; epithelial family extended beyond
## keratins/desmosomes to tight-junction/adhesion genes (Cldn/Ocln/Tjp/Epcam/Cdh1)
## since the batch1-only 8-gene signature surfaced two claudins (Cldn1, Cldn19).
BLOOD_REGEX <- "^(Hba|Hbb|Hbq|Hbz|Alas2|Slc4a1|Gypa|Rhd|Rhag|Ank1|Spta1|Sptb|Epb42)"
EPI_REGEX   <- "^(Krt|Dsg|Dsc|Pkp|Perp|Sprr|Flg|Lor|Ivl|Calml3|Slurp1|Lce|Cldn|Ocln|Tjp|Epcam|Cdh1)"

## ---- helpers (same as run_day7_day18_signatures.R) ----
strip_ensembl_version <- function(x) sub("\\.\\d+$", "", x)

log2cpm_tmm <- function(counts_mat) {
  dge <- DGEList(counts = round(pmax(counts_mat, 0)))
  dge <- calcNormFactors(dge, method = "TMM")
  cpm_mat <- cpm(dge)
  log2(cpm_mat + 1)
}

cpm_filter <- function(counts_mat, min_cpm = 1, min_frac = 0.25) {
  lib <- colSums(pmax(counts_mat, 0))
  lib[lib == 0] <- 1
  cpm_raw <- t(t(counts_mat) / lib) * 1e6
  min_samples <- max(2L, floor(ncol(counts_mat) * min_frac))
  keep <- rowSums(cpm_raw >= min_cpm) >= min_samples
  counts_mat[keep, , drop = FALSE]
}

zscore_rows <- function(mat) {
  zs <- t(scale(t(mat)))
  zs[is.na(zs)] <- 0
  zs
}

limma_tstat <- function(mat, y_factor, pos_label = "NonSurvivor") {
  y <- droplevels(factor(y_factor))
  grp_num <- as.integer(y == pos_label)
  design <- model.matrix(~ grp_num)
  fit <- suppressMessages(lmFit(mat, design))
  fit <- suppressMessages(eBayes(fit))
  tt <- fit$t[, 2]
  names(tt) <- rownames(mat)
  tt
}

make_loocv_folds <- function(y) {
  lapply(seq_along(y), function(i) list(train = setdiff(seq_along(y), i), test = i))
}

make_inner_foldid <- function(y, nfolds = 5L) {
  y <- factor(y)
  nfolds <- max(2L, min(as.integer(nfolds), min(table(y))))
  foldid <- integer(length(y))
  for (cl in levels(y)) {
    idx <- which(y == cl)
    foldid[idx] <- rep(seq_len(nfolds), length.out = length(idx))
  }
  foldid
}

make_repeated_stratified_folds <- function(y, K = 5L, repeats = 1L, seed = 123L) {
  y <- factor(y)
  cls <- levels(y)
  K_eff <- max(2L, min(as.integer(K), min(table(y))))
  set.seed(seed)
  folds <- list()
  f <- 1L
  for (r in seq_len(repeats)) {
    class_chunks <- lapply(cls, function(cl) {
      idx <- sample(which(y == cl))
      split(idx, cut(seq_along(idx), breaks = K_eff, labels = FALSE))
    })
    for (k in seq_len(K_eff)) {
      test_idx <- unlist(lapply(class_chunks, function(chlist) chlist[[k]]))
      folds[[f]] <- list(train = setdiff(seq_along(y), test_idx), test = test_idx)
      f <- f + 1L
    }
  }
  folds
}

map_ensembl_to_symbol <- function(ens_ids) {
  ens_ids <- sub("\\.\\d+$", "", ens_ids)
  is_mouse <- any(grepl("^ENSMUSG", ens_ids))
  sym <- rep(NA_character_, length(ens_ids))
  if (is_mouse && requireNamespace("org.Mm.eg.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    odb <- getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
    tab <- AnnotationDbi::select(odb, keys = ens_ids, keytype = "ENSEMBL", columns = "SYMBOL")
    sym[match(tab$ENSEMBL, ens_ids)] <- tab$SYMBOL
  }
  sym
}

normalize_for_contrast <- function(counts_mat, meta, do_batch_remove = TRUE) {
  counts_mat <- cpm_filter(counts_mat, min_cpm = 1, min_frac = 0.25)
  mat <- log2cpm_tmm(counts_mat)
  batch_removed <- FALSE
  if (do_batch_remove && "Batch" %in% names(meta) && nlevels(droplevels(factor(meta$Batch))) > 1L) {
    mat <- removeBatchEffect(mat, batch = droplevels(factor(meta$Batch)))
    batch_removed <- TRUE
  }
  mat <- mat[rowVars(mat) > 0, , drop = FALSE]
  list(mat = mat, mat_z = zscore_rows(mat), batch_removed = batch_removed)
}

fit_signature <- function(mat_z, y, top_genes, alpha) {
  M <- t(mat_z[top_genes, , drop = FALSE])
  foldid <- make_inner_foldid(y, nfolds = min(length(y), 5L))
  cv_fit <- suppressWarnings(cv.glmnet(
    M, y, alpha = alpha, family = "binomial",
    standardize = FALSE, foldid = foldid
  ))
  fit <- glmnet(
    M, y, alpha = alpha, lambda = cv_fit$lambda.min,
    family = "binomial", standardize = FALSE
  )
  coefs <- as.matrix(coef(fit))[-1, , drop = TRUE]
  nz <- which(coefs != 0)
  list(genes = names(coefs)[nz], coefs = coefs[nz], signature_size = length(nz))
}

predict_fold_alphas <- function(mat_z, y, tr, te, top_genes, alphas) {
  Mtr <- t(mat_z[top_genes, tr, drop = FALSE])
  Mte <- t(mat_z[top_genes, te, drop = FALSE])
  foldid <- make_inner_foldid(y[tr], nfolds = min(length(y[tr]), 5L))
  preds <- matrix(NA_real_, length(te), length(alphas))
  for (ai in seq_along(alphas)) {
    cv_fit <- suppressWarnings(cv.glmnet(
      Mtr, y[tr], alpha = alphas[ai], family = "binomial",
      standardize = FALSE, foldid = foldid
    ))
    fit <- glmnet(
      Mtr, y[tr], alpha = alphas[ai], lambda = cv_fit$lambda.min,
      family = "binomial", standardize = FALSE
    )
    preds[, ai] <- as.numeric(predict(fit, newx = Mte, type = "response"))
  }
  preds
}

loocv_predictions_by_alpha <- function(mat, mat_z, y, group_bin, K, alphas) {
  folds <- make_loocv_folds(y)
  preds <- matrix(NA_real_, length(y), length(alphas))
  for (fi in seq_along(folds)) {
    tr <- folds[[fi]]$train
    te <- folds[[fi]]$test
    t_tr <- limma_tstat(mat[, tr, drop = FALSE], group_bin[tr], pos_label = "NonSurvivor")
    top_tr <- names(t_tr)[order(abs(t_tr), decreasing = TRUE)[seq_len(min(K, length(t_tr)))]]
    preds[te, ] <- predict_fold_alphas(mat_z, y, tr, te, top_tr, alphas)
  }
  preds
}

rkfold_predictions_by_alpha <- function(mat, mat_z, y, group_bin, K, alphas,
                                        cv_k = RKFOLD_K, cv_repeats = RKFOLD_REPEATS) {
  pred_sum <- matrix(0, length(y), length(alphas))
  pred_cnt <- matrix(0L, length(y), length(alphas))
  for (rep_i in seq_len(cv_repeats)) {
    folds <- make_repeated_stratified_folds(y, K = cv_k, repeats = 1L, seed = 1000L + rep_i)
    for (fi in seq_along(folds)) {
      tr <- folds[[fi]]$train
      te <- folds[[fi]]$test
      t_tr <- limma_tstat(mat[, tr, drop = FALSE], group_bin[tr], pos_label = "NonSurvivor")
      top_tr <- names(t_tr)[order(abs(t_tr), decreasing = TRUE)[seq_len(min(K, length(t_tr)))]]
      fold_preds <- predict_fold_alphas(mat_z, y, tr, te, top_tr, alphas)
      for (j in seq_along(te)) {
        pred_sum[te[j], ] <- pred_sum[te[j], ] + fold_preds[j, ]
        pred_cnt[te[j], ] <- pred_cnt[te[j], ] + 1L
      }
    }
  }
  pred_sum / pmax(pred_cnt, 1L)
}

auc_from_preds <- function(y, preds) {
  as.numeric(auc(suppressWarnings(roc(y, preds, levels = c(0, 1), direction = "<", quiet = TRUE))))
}

roc_table <- function(y, preds) {
  roc_obj <- suppressWarnings(roc(y, preds, levels = c(0, 1), direction = "<", quiet = TRUE))
  rc <- as.data.frame(coords(roc_obj, "all",
    ret = c("threshold", "specificity", "sensitivity"), transpose = FALSE))
  rc$FPR <- 1 - rc$specificity
  rc$TPR <- rc$sensitivity
  list(table = rc, auc = as.numeric(auc(roc_obj)))
}

stability_flag <- function(loocv_auc, rkfold_auc) {
  if (is.na(loocv_auc) || is.na(rkfold_auc)) return("unknown")
  if (loocv_auc < 0.70 && rkfold_auc < 0.70) return("weak")
  if (loocv_auc >= 0.75 && abs(loocv_auc - rkfold_auc) <= 0.10) return("strong")
  if (loocv_auc >= 0.75 && rkfold_auc < loocv_auc - 0.10) return("loocv_only")
  "moderate"
}

## ----------------------------------------
## Load & preprocess counts (same as run_day7_day18_signatures.R)
## ----------------------------------------
meta_all <- read_tsv(INPUT_META, show_col_types = FALSE)
cts_tbl  <- read_tsv(INPUT_COUNTS, show_col_types = FALSE)
gene_col  <- names(cts_tbl)[1]
sample_cols <- setdiff(names(cts_tbl), gene_col)

ens_novers <- strip_ensembl_version(as.character(cts_tbl[[gene_col]]))
mat_raw <- as.matrix(cts_tbl[, sample_cols, drop = FALSE])
storage.mode(mat_raw) <- "numeric"
rownames(mat_raw) <- ens_novers
cts_df <- rowsum(mat_raw, group = rownames(mat_raw), reorder = FALSE)

ens_ids <- rownames(cts_df)
sym_vec <- map_ensembl_to_symbol(ens_ids)

keep <- !is.na(sym_vec) & nzchar(sym_vec)
cts_df_kept <- cts_df[keep, , drop = FALSE]
sym_kept <- sym_vec[keep]
mask <- !grepl("Rik$", sym_kept) & !grepl("^(Gm|Ps)", sym_kept)
cts_df_kept <- cts_df_kept[mask, , drop = FALSE]
sym_kept <- sym_kept[mask]
cts_sym <- rowsum(cts_df_kept, group = sym_kept, reorder = FALSE)

meta18 <- meta_all %>%
  filter(Timepoint == "Day18") %>%
  filter(Survival_Status %in% c("survivor", "non-survivor")) %>%
  mutate(
    Batch    = as.factor(as.character(Batch)),
    Treatment = tolower(Treatment),
    GroupBin  = if_else(tolower(Survival_Status) == "non-survivor", "NonSurvivor", "Survivor")
  )
meta18_untreated_batch1 <- meta18 %>%
  filter(Treatment == "n", as.character(Batch) == TARGET_BATCH) %>%
  transmute(Core_Sample_Id, Batch, Treatment, GroupBin)

message(sprintf("Day18_Untreated_SvNS restricted to %s: N=%d (%s)",
                TARGET_BATCH, nrow(meta18_untreated_batch1),
                paste(names(table(meta18_untreated_batch1$GroupBin)),
                      table(meta18_untreated_batch1$GroupBin), sep = "=", collapse = ", ")))

## ==========================================================
## Day18_Untreated_SvNS, batch_1 only, contamination excluded, ~15 genes
## ==========================================================
tag <- "Day18_Untreated_SvNS"
message("\n=== ", tag, " (batch_1-only, contam-filtered, ~15-gene target) ===")

common_ids <- intersect(colnames(cts_sym), meta18_untreated_batch1$Core_Sample_Id)
Xraw <- cts_sym[, common_ids, drop = FALSE]
meta <- as.data.frame(meta18_untreated_batch1[match(common_ids, meta18_untreated_batch1$Core_Sample_Id), , drop = FALSE])
y <- as.integer(meta$GroupBin == "NonSurvivor")

norm  <- normalize_for_contrast(Xraw, meta, do_batch_remove = FALSE)  # single batch present; no removal needed
mat   <- norm$mat
mat_z <- norm$mat_z
message(sprintf("  Genes after CPM filter: %d; batch_removed=%s; N=%d",
                nrow(mat), norm$batch_removed, length(y)))

## ---- exclude contamination gene families before any ranking ----
is_blood <- grepl(BLOOD_REGEX, rownames(mat))
is_epi   <- grepl(EPI_REGEX, rownames(mat))
is_contam <- is_blood | is_epi

excl_df <- data.frame(
  gene = rownames(mat)[is_contam],
  family = ifelse(is_blood[is_contam], "blood_erythroid", "epithelial_tightjunction")
)
excl_df <- excl_df[order(excl_df$family, excl_df$gene), ]
write.csv(excl_df, file.path(OUTDIR, "excluded_contamination_genes.csv"), row.names = FALSE)
message(sprintf("  Excluding %d contamination-family genes (%d blood/erythroid, %d epithelial/tight-junction)",
                nrow(excl_df), sum(is_blood), sum(is_epi)))

mat   <- mat[!is_contam, , drop = FALSE]
mat_z <- mat_z[!is_contam, , drop = FALSE]

## ---- K/alpha sweep on the filtered gene universe, 12-18 gene target band ----
full_grid <- data.frame()
for (K in K_GRID) {
  loocv_preds  <- loocv_predictions_by_alpha(mat, mat_z, y, meta$GroupBin, K, ALPHAS)
  rkfold_preds <- rkfold_predictions_by_alpha(mat, mat_z, y, meta$GroupBin, K, ALPHAS)

  loocv_auc  <- apply(loocv_preds,  2, auc_from_preds, y = y)
  rkfold_auc <- apply(rkfold_preds, 2, auc_from_preds, y = y)

  tstat_all <- limma_tstat(mat, meta$GroupBin, pos_label = "NonSurvivor")
  top_all <- names(tstat_all)[order(abs(tstat_all), decreasing = TRUE)[seq_len(min(K, length(tstat_all)))]]
  sig_sizes <- vapply(ALPHAS, function(alpha) {
    fit_signature(mat_z, y, top_all, alpha)$signature_size
  }, integer(1))

  k_df <- data.frame(
    K = K, alpha = ALPHAS,
    AUC_LOOCV = loocv_auc, AUC_repeated_kfold = rkfold_auc,
    SigSize = sig_sizes,
    Within_constraint = sig_sizes >= SIG_MIN & sig_sizes <= SIG_MAX
  )
  write.csv(k_df, file.path(OUTDIR, sprintf("K%03d_AUC_by_alpha.csv", K)), row.names = FALSE)
  full_grid <- rbind(full_grid, k_df)

  message(sprintf(
    "    K=%d done. Best LOOCV alpha=%.2f AUC=%.4f; best rkfold alpha=%.2f AUC=%.4f",
    K, ALPHAS[which.max(loocv_auc)], max(loocv_auc, na.rm = TRUE),
    ALPHAS[which.max(rkfold_auc)], max(rkfold_auc, na.rm = TRUE)
  ))
}

full_grid <- full_grid[order(-full_grid$AUC_LOOCV, -full_grid$AUC_repeated_kfold, full_grid$SigSize), ]
write.csv(full_grid, file.path(OUTDIR, "full_auc_grid.csv"), row.names = FALSE)

candidates <- full_grid[full_grid$Within_constraint, , drop = FALSE]
if (nrow(candidates) == 0) {
  stop(sprintf("No %s candidates satisfy the %d-%d gene constraint after contamination filtering.",
               tag, SIG_MIN, SIG_MAX))
}
best_row <- candidates[order(-candidates$AUC_LOOCV, -candidates$AUC_repeated_kfold, candidates$SigSize)[1], ]
best_unconstrained_row <- full_grid[order(-full_grid$AUC_LOOCV, -full_grid$AUC_repeated_kfold, full_grid$SigSize)[1], ]

best_K     <- best_row$K
best_alpha <- best_row$alpha

selected_loocv  <- loocv_predictions_by_alpha(mat, mat_z, y, meta$GroupBin, best_K, best_alpha)[, 1]
selected_rkfold <- rkfold_predictions_by_alpha(mat, mat_z, y, meta$GroupBin, best_K, best_alpha)[, 1]
roc_loocv <- roc_table(y, selected_loocv)

tstat_best <- limma_tstat(mat, meta$GroupBin, pos_label = "NonSurvivor")
top_best <- names(tstat_best)[order(abs(tstat_best), decreasing = TRUE)[seq_len(min(best_K, length(tstat_best)))]]
sig_best <- fit_signature(mat_z, y, top_best, best_alpha)

signature_tbl <- data.frame(
  gene = sig_best$genes, coef = as.numeric(sig_best$coefs),
  direction = ifelse(sig_best$coefs >= 0, "NS_high", "S_high"),
  K = best_K, alpha = best_alpha,
  AUC_LOOCV = roc_loocv$auc,
  AUC_repeated_kfold = auc_from_preds(y, selected_rkfold),
  stringsAsFactors = FALSE
)
signature_tbl <- signature_tbl[order(abs(signature_tbl$coef), decreasing = TRUE), ]

write.csv(signature_tbl, file.path(OUTDIR, "selected_signature_genes.csv"), row.names = FALSE)
write.csv(roc_loocv$table, file.path(OUTDIR, "selected_ROC_data_loocv.csv"), row.names = FALSE)
write.csv(
  data.frame(Core_Sample_Id = meta$Core_Sample_Id, Group = meta$GroupBin, LOOCV_score = selected_loocv),
  file.path(OUTDIR, "selected_LOOCV_scores.csv"), row.names = FALSE
)
write.csv(
  data.frame(Core_Sample_Id = meta$Core_Sample_Id, Group = meta$GroupBin, rkfold_score = selected_rkfold),
  file.path(OUTDIR, "selected_rkfold_scores.csv"), row.names = FALSE
)

selected <- data.frame(
  Contrast = tag,
  normalization = "log2cpm_tmm_batch1_only_contam_filtered",
  Batch_removed = norm$batch_removed,
  N = length(y),
  Selected_K = best_K, Selected_alpha = best_alpha,
  AUC_LOOCV = roc_loocv$auc,
  AUC_repeated_kfold = auc_from_preds(y, selected_rkfold),
  Signature_genes = sig_best$signature_size,
  Selected_by_constraint = TRUE,
  Stability_flag = stability_flag(roc_loocv$auc, auc_from_preds(y, selected_rkfold)),
  Best_unconstrained_K = best_unconstrained_row$K,
  Best_unconstrained_alpha = best_unconstrained_row$alpha,
  Best_unconstrained_genes = best_unconstrained_row$SigSize,
  Best_unconstrained_within_constraint = best_unconstrained_row$Within_constraint,
  Best_LOOCV_unconstrained = best_unconstrained_row$AUC_LOOCV,
  Best_rkfold_unconstrained = best_unconstrained_row$AUC_repeated_kfold,
  stringsAsFactors = FALSE
)
write.csv(selected, file.path(OUTDIR, "selected_best_candidate.csv"), row.names = FALSE)

metrics <- c(
  sprintf("=== %s: batch_1-only, contam-filtered, ~15-gene target (N=%d) ===", tag, length(y)),
  sprintf("normalization=log2cpm_tmm (batch_1 only, blood+epithelial/tight-junction genes excluded pre-ranking), K=%d, alpha=%.2f",
          best_K, best_alpha),
  sprintf("LOOCV AUC: %.4f", selected$AUC_LOOCV),
  sprintf("Repeated %dx%d-fold AUC: %.4f", RKFOLD_REPEATS, RKFOLD_K, selected$AUC_repeated_kfold),
  sprintf("Signature size: %d genes", selected$Signature_genes),
  sprintf("Stability flag: %s", selected$Stability_flag),
  sprintf("Contamination genes excluded from candidacy: %d", nrow(excl_df)),
  "",
  "Same pipeline as run_day7_day18_signatures.R (CPM pre-filter + TMM log2CPM + per-fold limma ranking), applied to",
  "Day18_Untreated_SvNS restricted to batch_1 only (N=16), with blood/erythroid and",
  "epithelial/tight-junction gene families (Krt/Dsg/Dsc/Pkp/Perp/... + Cldn/Ocln/Tjp/Epcam/Cdh1)",
  "removed from the candidate pool before ranking at every K, and the signature-size selection",
  "constraint narrowed from the default (8-30) to 12-18 genes to target a larger, ~15-gene panel."
)
writeLines(metrics, file.path(OUTDIR, "metrics_summary.txt"))
cat(paste0(metrics, collapse = "\n"), "\n")

## ==========================================================
## DeLong CI + permutation test (scoped to this one contrast)
## ==========================================================
message(sprintf("\n  Running %d-permutation test (fixed K=%d, alpha=%.2f)...", N_PERM, best_K, best_alpha))

roc_obj <- suppressWarnings(roc(y, selected_loocv, levels = c(0, 1), direction = "<", quiet = TRUE))
ci_delong <- tryCatch(
  as.numeric(ci.auc(roc_obj, method = "delong", conf.level = 0.95)),
  error = function(e) c(NA, roc_loocv$auc, NA)
)

loocv_fixed_k_alpha <- function(mat, mat_z, y, group_bin, K, alpha) {
  n <- length(y)
  preds <- numeric(n)
  for (i in seq_len(n)) {
    tr <- setdiff(seq_len(n), i); te <- i
    t_tr <- limma_tstat(mat[, tr, drop = FALSE], group_bin[tr], "NonSurvivor")
    top_genes <- names(t_tr)[order(abs(t_tr), decreasing = TRUE)[seq_len(min(K, length(t_tr)))]]
    Mtr <- t(mat_z[top_genes, tr, drop = FALSE])
    Mte <- t(mat_z[top_genes, te, drop = FALSE])
    fid <- make_inner_foldid(y[tr], nfolds = min(length(tr), 5L))
    cv_fit <- suppressWarnings(cv.glmnet(Mtr, y[tr], alpha = alpha, family = "binomial",
                                         standardize = FALSE, foldid = fid))
    fit <- glmnet(Mtr, y[tr], alpha = alpha, lambda = cv_fit$lambda.min,
                  family = "binomial", standardize = FALSE)
    preds[te] <- as.numeric(predict(fit, newx = Mte, type = "response"))
  }
  preds
}

null_aucs <- numeric(N_PERM)
set.seed(42)
for (perm_i in seq_len(N_PERM)) {
  y_perm <- sample(y)
  gb_perm <- ifelse(y_perm == 1, "NonSurvivor", "Survivor")
  preds_perm <- tryCatch(
    loocv_fixed_k_alpha(mat, mat_z, y_perm, gb_perm, best_K, best_alpha),
    error = function(e) rep(0.5, length(y_perm))
  )
  null_aucs[perm_i] <- auc_from_preds(y_perm, preds_perm)
  if (perm_i %% 50 == 0) message(sprintf("    perm %d/%d done", perm_i, N_PERM))
}
perm_p <- mean(null_aucs >= roc_loocv$auc)
message(sprintf("  Permutation p-value: %.4f (obs AUC=%.4f, null median=%.4f)",
                perm_p, roc_loocv$auc, median(null_aucs)))

val_dir <- file.path(OUTDIR, "validation_ci_perm")
dir.create(val_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(perm = seq_len(N_PERM), null_auc = null_aucs),
          file.path(val_dir, "permutation_null_aucs.csv"), row.names = FALSE)

png(file.path(val_dir, "permutation_null_histogram.png"), width = 700, height = 500, res = 96)
hist(null_aucs, breaks = 30, col = "steelblue", border = "white",
     main = sprintf("%s (batch_1, contam-filtered, ~15-gene)\nNull AUC distribution (N=%d permutations)", tag, N_PERM),
     xlab = "Permuted AUC", xlim = c(0, 1))
abline(v = roc_loocv$auc, col = "red", lwd = 2.5)
legend("topleft", legend = sprintf("Observed AUC = %.3f\np = %.4f", roc_loocv$auc, perm_p),
       col = "red", lwd = 2, bty = "n")
dev.off()

write.csv(
  data.frame(Contrast = tag, N = length(y), LOOCV_AUC = roc_loocv$auc,
             DeLong_CI_lo = ci_delong[1], DeLong_CI_hi = ci_delong[3],
             Perm_p = perm_p, N_perms = N_PERM),
  file.path(val_dir, "validation_ci_perm_summary.csv"), row.names = FALSE
)

## ==========================================================
## Comparison vs. pooled and vs. batch1-only (unfiltered, 8-gene)
## ==========================================================
orig_best  <- read.csv(file.path(SIG_ROOT, tag, "selected_best_candidate.csv"))
orig_genes <- read.csv(file.path(SIG_ROOT, tag, "selected_signature_genes.csv"))$gene
orig_vsum  <- tryCatch(
  read.csv(file.path(SIG_ROOT, "validation_ci_perm_summary.csv")) %>% filter(Contrast == tag),
  error = function(e) data.frame(Perm_p = NA)
)
orig_n <- nrow(read.csv(file.path(SIG_ROOT, tag, "selected_LOOCV_scores.csv")))

b1_best  <- read.csv(file.path(BATCH1_ROOT, tag, "selected_best_candidate.csv"))
b1_genes <- read.csv(file.path(BATCH1_ROOT, tag, "selected_signature_genes.csv"))$gene
b1_vsum  <- read.csv(file.path(BATCH1_ROOT, tag, "validation_ci_perm", "validation_ci_perm_summary.csv"))

comparison <- data.frame(
  Version = c("pooled_batches", "batch1_only", "batch1_contam_filtered_15gene"),
  N = c(orig_n, b1_best$N[1], length(y)),
  K = c(orig_best$Selected_K[1], b1_best$Selected_K[1], best_K),
  alpha = c(orig_best$Selected_alpha[1], b1_best$Selected_alpha[1], best_alpha),
  N_genes = c(length(orig_genes), length(b1_genes), nrow(signature_tbl)),
  LOOCV_AUC = c(orig_best$AUC_LOOCV[1], b1_best$AUC_LOOCV[1], roc_loocv$auc),
  rkfold_AUC = c(orig_best$AUC_repeated_kfold[1], b1_best$AUC_repeated_kfold[1], selected$AUC_repeated_kfold),
  DeLong_lo = c(NA, b1_vsum$DeLong_CI_lo[1], ci_delong[1]),
  DeLong_hi = c(NA, b1_vsum$DeLong_CI_hi[1], ci_delong[3]),
  Perm_p = c(if (nrow(orig_vsum) > 0) orig_vsum$Perm_p[1] else NA, b1_vsum$Perm_p[1], perm_p),
  Genes = c(paste(orig_genes, collapse = ", "), paste(b1_genes, collapse = ", "),
            paste(signature_tbl$gene, collapse = ", "))
)
write.csv(comparison, file.path(OUTDIR, "comparison_vs_pooled_and_batch1.csv"), row.names = FALSE)

message("\n=== COMPARISON: pooled vs batch1-only (8-gene) vs batch1 contam-filtered ~15-gene ===")
print(comparison[, c("Version","N","K","alpha","N_genes","LOOCV_AUC","rkfold_AUC","DeLong_lo","DeLong_hi","Perm_p")])
message("\nNew signature genes: ", paste(signature_tbl$gene, collapse = ", "))
message("\nAll outputs saved to: ", normalizePath(OUTDIR))
