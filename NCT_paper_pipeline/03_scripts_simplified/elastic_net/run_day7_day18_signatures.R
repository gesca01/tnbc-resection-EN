## Elastic-net survival signatures for the six day-7 / day-18 contrasts
## (treated, untreated, combined).
##   - CPM filter (CPM >= 1 in >= 25% of samples), TMM, log2CPM
##   - limma batch removal when a contrast spans two sequencing batches
##   - row z-scoring
##   - limma moderated-t gene ranking inside each leave-one-out fold
##   - grid over K (number of ranked genes) and alpha, final signature 8-30 genes

options(stringsAsFactors = FALSE)
set.seed(123)

.libPaths(c("./NCT_paper_pipeline/renv_lib", .libPaths()))
source("./NCT_paper_pipeline/03_scripts_core/elastic_net/elastic_net_config.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(matrixStats)
  library(limma)
  library(glmnet)
  library(pROC)
  library(edgeR)
})

OUTDIR <- file.path(
  PROJECT_ROOT,
  "02_cache",
  "elastic_net_day7_day18_source_match",
  "results_EN_signatures"
)

K_GRID  <- c(5, 7, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200)
ALPHAS  <- c(0.01, 0.02, 0.03, seq(0.05, 1.0, by = 0.05))
SIG_MIN <- 8L
SIG_MAX <- 30L
RKFOLD_K     <- 5L
RKFOLD_REPEATS <- 10L

safe_dir <- function(...) {
  d <- file.path(...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

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
  is_human <- any(grepl("^ENSG", ens_ids))
  sym <- rep(NA_character_, length(ens_ids))

  if (is_mouse && requireNamespace("org.Mm.eg.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    odb <- getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
    tab <- AnnotationDbi::select(odb, keys = ens_ids, keytype = "ENSEMBL", columns = "SYMBOL")
    sym[match(tab$ENSEMBL, ens_ids)] <- tab$SYMBOL
  } else if (is_human && requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
             requireNamespace("AnnotationDbi", quietly = TRUE)) {
    odb <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
    tab <- AnnotationDbi::select(odb, keys = ens_ids, keytype = "ENSEMBL", columns = "SYMBOL")
    sym[match(tab$ENSEMBL, ens_ids)] <- tab$SYMBOL
  } else if (requireNamespace("biomaRt", quietly = TRUE)) {
    host_url <- "https://dec2021.archive.ensembl.org"
    if (is_mouse) {
      mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host_url)
      tab <- biomaRt::getBM(
        attributes = c("ensembl_gene_id", "mgi_symbol"),
        filters = "ensembl_gene_id", values = ens_ids, mart = mart
      )
      sym[match(tab$ensembl_gene_id, ens_ids)] <- tab$mgi_symbol
    } else if (is_human) {
      mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host_url)
      tab <- biomaRt::getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol"),
        filters = "ensembl_gene_id", values = ens_ids, mart = mart
      )
      sym[match(tab$ensembl_gene_id, ens_ids)] <- tab$hgnc_symbol
    }
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

run_binary <- function(tag, counts_sym, meta_sub, do_batch_remove = TRUE,
                           k_grid = K_GRID, alphas = ALPHAS,
                           sig_min = SIG_MIN, sig_max = SIG_MAX) {
  set.seed(123)
  message("\n=== ", tag, " ===")

  common_ids <- intersect(colnames(counts_sym), meta_sub$Core_Sample_Id)
  stopifnot(length(common_ids) >= 2)
  Xraw <- counts_sym[, common_ids, drop = FALSE]
  meta <- as.data.frame(meta_sub[match(common_ids, meta_sub$Core_Sample_Id), , drop = FALSE])
  y <- as.integer(meta$GroupBin == "NonSurvivor")
  outdir <- safe_dir(OUTDIR, tag)

  norm <- normalize_for_contrast(Xraw, meta, do_batch_remove)
  mat   <- norm$mat
  mat_z <- norm$mat_z
  message(sprintf("  Genes after CPM filter: %d; batch_removed=%s", nrow(mat), norm$batch_removed))

  full_grid <- data.frame()

  for (K in k_grid) {
    loocv_preds  <- loocv_predictions_by_alpha(mat, mat_z, y, meta$GroupBin, K, alphas)
    rkfold_preds <- rkfold_predictions_by_alpha(mat, mat_z, y, meta$GroupBin, K, alphas)

    loocv_auc  <- apply(loocv_preds,  2, auc_from_preds, y = y)
    rkfold_auc <- apply(rkfold_preds, 2, auc_from_preds, y = y)

    tstat_all <- limma_tstat(mat, meta$GroupBin, pos_label = "NonSurvivor")
    top_all <- names(tstat_all)[order(abs(tstat_all), decreasing = TRUE)[seq_len(min(K, length(tstat_all)))]]
    sig_sizes <- vapply(alphas, function(alpha) {
      fit_signature(mat_z, y, top_all, alpha)$signature_size
    }, integer(1))

    k_df <- data.frame(
      K = K, alpha = alphas,
      AUC_LOOCV = loocv_auc, AUC_repeated_kfold = rkfold_auc,
      SigSize = sig_sizes,
      Within_constraint = sig_sizes >= sig_min & sig_sizes <= sig_max
    )
    write.csv(k_df, file.path(outdir, sprintf("K%03d_AUC_by_alpha.csv", K)), row.names = FALSE)
    full_grid <- rbind(full_grid, k_df)

    message(sprintf(
      "    K=%d done. Best LOOCV alpha=%.2f AUC=%.4f; best rkfold alpha=%.2f AUC=%.4f",
      K, alphas[which.max(loocv_auc)], max(loocv_auc, na.rm = TRUE),
      alphas[which.max(rkfold_auc)], max(rkfold_auc, na.rm = TRUE)
    ))
  }

  full_grid <- full_grid[order(-full_grid$AUC_LOOCV, -full_grid$AUC_repeated_kfold, full_grid$SigSize), ]
  write.csv(full_grid, file.path(outdir, "full_auc_grid.csv"), row.names = FALSE)

  candidates <- full_grid[full_grid$Within_constraint, , drop = FALSE]
  if (nrow(candidates) == 0) {
    stop(sprintf("No %s candidates satisfy the %d-%d gene constraint.", tag, sig_min, sig_max))
  } else {
    best_row <- candidates[order(-candidates$AUC_LOOCV, -candidates$AUC_repeated_kfold, candidates$SigSize)[1], ]
    selected_by_constraint <- TRUE
  }
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

  write.csv(signature_tbl, file.path(outdir, "selected_signature_genes.csv"), row.names = FALSE)
  write.csv(roc_loocv$table, file.path(outdir, "selected_ROC_data_loocv.csv"), row.names = FALSE)
  write.csv(
    data.frame(Core_Sample_Id = meta$Core_Sample_Id, Group = meta$GroupBin, LOOCV_score = selected_loocv),
    file.path(outdir, "selected_LOOCV_scores.csv"), row.names = FALSE
  )
  write.csv(
    data.frame(Core_Sample_Id = meta$Core_Sample_Id, Group = meta$GroupBin, rkfold_score = selected_rkfold),
    file.path(outdir, "selected_rkfold_scores.csv"), row.names = FALSE
  )

  selected <- data.frame(
    Contrast = tag,
    normalization = "log2cpm_tmm",
    Batch_removed = norm$batch_removed,
    Selected_K = best_K, Selected_alpha = best_alpha,
    AUC_LOOCV = roc_loocv$auc,
    AUC_repeated_kfold = auc_from_preds(y, selected_rkfold),
    Signature_genes = sig_best$signature_size,
    Selected_by_constraint = selected_by_constraint,
    Stability_flag = stability_flag(roc_loocv$auc, auc_from_preds(y, selected_rkfold)),
    Best_unconstrained_K = best_unconstrained_row$K,
    Best_unconstrained_alpha = best_unconstrained_row$alpha,
    Best_unconstrained_genes = best_unconstrained_row$SigSize,
    Best_unconstrained_within_constraint = best_unconstrained_row$Within_constraint,
    Best_LOOCV_unconstrained = best_unconstrained_row$AUC_LOOCV,
    Best_rkfold_unconstrained = best_unconstrained_row$AUC_repeated_kfold,
    stringsAsFactors = FALSE
  )
  write.csv(selected, file.path(outdir, "selected_best_candidate.csv"), row.names = FALSE)

  metrics <- c(
    sprintf("=== %s ===", tag),
    sprintf("normalization=log2cpm_tmm, K=%d, alpha=%.2f", best_K, best_alpha),
    sprintf("LOOCV AUC: %.4f", selected$AUC_LOOCV),
    sprintf("Repeated %dx%d-fold AUC: %.4f", RKFOLD_REPEATS, RKFOLD_K, selected$AUC_repeated_kfold),
    sprintf("Signature size: %d genes", selected$Signature_genes),
    sprintf("Stability flag: %s", selected$Stability_flag),
    "",
    "TMM normalization, CPM pre-filter (CPM>=1 in >=25% samples), 8-30 gene signature.",
    "Selection: highest LOOCV AUC within the size limit; ties broken by repeated k-fold AUC, then smaller signature."
  )
  writeLines(metrics, file.path(outdir, "metrics_summary.txt"))
  cat(paste0(metrics, collapse = "\n"), "\n")

  selected
}

write_summary <- function(summary_df) {
  summary_df$LOOCV_minus_rkfold <- summary_df$AUC_LOOCV - summary_df$AUC_repeated_kfold
  write.csv(summary_df, file.path(OUTDIR, "signature_summary.csv"), row.names = FALSE)

  fmt <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))
  rows <- apply(summary_df, 1, function(row) {
    sprintf("| %s | %s | %s | %s | %s | %s |",
      row[["Contrast"]],
      fmt(as.numeric(row[["AUC_LOOCV"]])),
      fmt(as.numeric(row[["AUC_repeated_kfold"]])),
      row[["Selected_K"]], sprintf("%.2f", as.numeric(row[["Selected_alpha"]])),
      row[["Signature_genes"]])
  })
  md <- c(
    "# Day 7 / day 18 elastic-net signatures",
    "",
    "| Contrast | LOOCV AUC | Repeated 10x5-fold AUC | K | alpha | Genes |",
    "|---|---:|---:|---:|---:|---:|",
    rows
  )
  writeLines(md, file.path(OUTDIR, "signature_summary.md"))
  invisible(summary_df)
}

## ----------------------------------------
## Load & preprocess counts
## ----------------------------------------
safe_dir(OUTDIR)

meta_all <- read_tsv(INPUT_METADATA, show_col_types = FALSE)
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
is_mouse <- any(grepl("^ENSMUSG", ens_ids))
is_human <- any(grepl("^ENSG", ens_ids))

if (all(is.na(sym_vec))) {
  message("WARNING: Symbol mapping failed. Proceeding with Ensembl IDs.")
  cts_sym <- cts_df
} else {
  keep <- !is.na(sym_vec) & nzchar(sym_vec)
  cts_df_kept <- cts_df[keep, , drop = FALSE]
  sym_kept <- sym_vec[keep]
  if (is_mouse) {
    mask <- !grepl("Rik$", sym_kept) & !grepl("^(Gm|Ps)", sym_kept)
    cts_df_kept <- cts_df_kept[mask, , drop = FALSE]
    sym_kept <- sym_kept[mask]
  } else if (is_human) {
    mask <- !grepl("RIK$", sym_kept) & !grepl("^(GM|PS)", sym_kept)
    cts_df_kept <- cts_df_kept[mask, , drop = FALSE]
    sym_kept <- sym_kept[mask]
  }
  cts_sym <- rowsum(cts_df_kept, group = sym_kept, reorder = FALSE)
}

mk_meta_day <- function(day) {
  meta_all %>%
    filter(Timepoint == day) %>%
    filter(Survival_Status %in% c("survivor", "non-survivor")) %>%
    mutate(
      Batch    = as.factor(as.character(Batch)),
      Treatment = tolower(Treatment),
      GroupBin  = if_else(tolower(Survival_Status) == "non-survivor", "NonSurvivor", "Survivor")
    )
}

meta18 <- mk_meta_day("Day18")
meta07 <- mk_meta_day("Day7")

meta18_combined   <- meta18 %>% transmute(Core_Sample_Id, Batch, Treatment, GroupBin)
meta18_treated    <- meta18 %>% filter(Treatment == "y") %>% transmute(Core_Sample_Id, Batch, Treatment, GroupBin)
meta18_untreated  <- meta18 %>% filter(Treatment == "n") %>% transmute(Core_Sample_Id, Batch, Treatment, GroupBin)
meta07_combined   <- meta07 %>% transmute(Core_Sample_Id, Batch, Treatment, GroupBin)
meta07_treated    <- meta07 %>% filter(Treatment == "y") %>% transmute(Core_Sample_Id, Batch, Treatment, GroupBin)
meta07_untreated  <- meta07 %>% filter(Treatment == "n") %>% transmute(Core_Sample_Id, Batch, Treatment, GroupBin)

summary_df <- bind_rows(
  run_binary("Day18_Combined_SvNS",   cts_sym, meta18_combined,  do_batch_remove = TRUE),
  run_binary("Day18_Treated_SvNS",    cts_sym, meta18_treated,   do_batch_remove = TRUE),
  run_binary("Day18_Untreated_SvNS",  cts_sym, meta18_untreated, do_batch_remove = TRUE),
  run_binary("Day7_Combined_SvNS",    cts_sym, meta07_combined,  do_batch_remove = TRUE),
  run_binary("Day7_Treated_SvNS",     cts_sym, meta07_treated,   do_batch_remove = TRUE),
  run_binary("Day7_Untreated_SvNS",   cts_sym, meta07_untreated, do_batch_remove = TRUE)
)

summary_df <- write_summary(summary_df)

cat("\nOutputs saved to:\n")
cat("  ", normalizePath(OUTDIR), "\n")
