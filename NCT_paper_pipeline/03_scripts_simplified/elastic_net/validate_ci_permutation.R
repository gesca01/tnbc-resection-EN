## DeLong 95% CI and label-permutation test for the day-7 / day-18 signatures.
## Each permutation shuffles outcome labels and reruns the full per-fold pipeline
## at the selected K and alpha. Number of permutations: N_PERM (default 200, as in the paper).

options(stringsAsFactors = FALSE)
set.seed(42)

.libPaths(c(
  "./NCT_paper_pipeline/renv_lib",
  .libPaths()
))

PROJECT_ROOT <- "./NCT_paper_pipeline"
source(file.path(PROJECT_ROOT, "03_scripts_core", "elastic_net", "elastic_net_config.R"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(matrixStats)
  library(limma)
  library(glmnet)
  library(pROC)
  library(edgeR)
})

OUTDIR <- file.path(PROJECT_ROOT, "02_cache",
  "elastic_net_day7_day18_source_match", "results_EN_signatures")
N_PERM <- as.integer(Sys.getenv("N_PERM", "200"))
CONTRASTS <- c("Day18_Combined_SvNS", "Day18_Treated_SvNS", "Day18_Untreated_SvNS",
               "Day7_Combined_SvNS", "Day7_Treated_SvNS", "Day7_Untreated_SvNS")

## ---- Helpers (same as run_day7_day18_signatures.R) ----------------------------------
strip_ensembl_version <- function(x) sub("\\.\\d+$", "", x)

log2cpm_tmm <- function(counts_mat) {
  dge <- DGEList(counts = round(pmax(counts_mat, 0)))
  dge <- calcNormFactors(dge, method = "TMM")
  log2(cpm(dge) + 1)
}

cpm_filter <- function(counts_mat, min_cpm = 1, min_frac = 0.25) {
  lib <- colSums(pmax(counts_mat, 0))
  lib[lib == 0] <- 1
  cpm_raw <- t(t(counts_mat) / lib) * 1e6
  min_samples <- max(2L, floor(ncol(counts_mat) * min_frac))
  counts_mat[rowSums(cpm_raw >= min_cpm) >= min_samples, , drop = FALSE]
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
  fit <- suppressMessages(eBayes(lmFit(mat, design)))
  setNames(fit$t[, 2], rownames(mat))
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

normalize_contrast <- function(counts_mat, meta, do_batch_remove = TRUE) {
  counts_mat <- cpm_filter(counts_mat, min_cpm = 1, min_frac = 0.25)
  mat <- log2cpm_tmm(counts_mat)
  if (do_batch_remove && "Batch" %in% names(meta) &&
      nlevels(droplevels(factor(meta$Batch))) > 1L) {
    mat <- removeBatchEffect(mat, batch = droplevels(factor(meta$Batch)))
  }
  mat <- mat[rowVars(mat) > 0, , drop = FALSE]
  list(mat = mat, mat_z = zscore_rows(mat))
}

loocv_fixed_k_alpha <- function(mat, mat_z, y, group_bin, K, alpha) {
  n <- length(y)
  preds <- numeric(n)
  for (i in seq_len(n)) {
    tr <- setdiff(seq_len(n), i)
    te <- i
    t_tr <- limma_tstat(mat[, tr, drop = FALSE], group_bin[tr], "NonSurvivor")
    top_genes <- names(t_tr)[order(abs(t_tr), decreasing = TRUE)[seq_len(min(K, length(t_tr)))]]
    Mtr <- t(mat_z[top_genes, tr, drop = FALSE])
    Mte <- t(mat_z[top_genes, te, drop = FALSE])
    fid <- make_inner_foldid(y[tr], nfolds = min(length(tr), 5L))
    cv_fit <- suppressWarnings(cv.glmnet(
      Mtr, y[tr], alpha = alpha, family = "binomial",
      standardize = FALSE, foldid = fid
    ))
    fit <- glmnet(Mtr, y[tr], alpha = alpha, lambda = cv_fit$lambda.min,
                  family = "binomial", standardize = FALSE)
    preds[te] <- as.numeric(predict(fit, newx = Mte, type = "response"))
  }
  preds
}

auc_roc <- function(y, preds) {
  as.numeric(auc(suppressWarnings(
    roc(y, preds, levels = c(0, 1), direction = "<", quiet = TRUE)
  )))
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

## ---- Load raw data ------------------------------------------
message("Loading raw counts and metadata...")
INPUT_COUNTS   <- file.path(PROJECT_ROOT, "01_inputs", "merged_counts.tsv")
INPUT_METADATA <- file.path(PROJECT_ROOT, "01_inputs", "merged_metadata_all.tsv")

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
if (all(is.na(sym_vec))) {
  message("WARNING: Symbol mapping failed. Using Ensembl IDs.")
  cts_sym <- cts_df
} else {
  keep <- !is.na(sym_vec) & nzchar(sym_vec)
  cts_df_k <- cts_df[keep, , drop = FALSE]
  sym_k <- sym_vec[keep]
  is_mouse <- any(grepl("^ENSMUSG", ens_ids))
  if (is_mouse) {
    mask <- !grepl("Rik$", sym_k) & !grepl("^(Gm|Ps)", sym_k)
    cts_df_k <- cts_df_k[mask, , drop = FALSE]
    sym_k <- sym_k[mask]
  }
  cts_sym <- rowsum(cts_df_k, group = sym_k, reorder = FALSE)
}

mk_meta_day <- function(day) {
  meta_all %>%
    filter(Timepoint == day) %>%
    filter(Survival_Status %in% c("survivor", "non-survivor")) %>%
    mutate(
      Batch     = as.factor(as.character(Batch)),
      Treatment = tolower(Treatment),
      GroupBin  = if_else(tolower(Survival_Status) == "non-survivor", "NonSurvivor", "Survivor")
    )
}

meta18 <- mk_meta_day("Day18")
meta07 <- mk_meta_day("Day7")

contrast_meta <- list(
  Day18_Combined_SvNS  = meta18 %>% transmute(Core_Sample_Id, Batch, GroupBin),
  Day18_Treated_SvNS   = meta18 %>% filter(Treatment == "y") %>% transmute(Core_Sample_Id, Batch, GroupBin),
  Day18_Untreated_SvNS = meta18 %>% filter(Treatment == "n") %>% transmute(Core_Sample_Id, Batch, GroupBin),
  Day7_Combined_SvNS   = meta07 %>% transmute(Core_Sample_Id, Batch, GroupBin),
  Day7_Treated_SvNS    = meta07 %>% filter(Treatment == "y") %>% transmute(Core_Sample_Id, Batch, GroupBin),
  Day7_Untreated_SvNS  = meta07 %>% filter(Treatment == "n") %>% transmute(Core_Sample_Id, Batch, GroupBin)
)

## ---- Run validation per contrast ----------------------------
results_all <- list()

for (contrast in CONTRASTS) {
  message("\n=== ", contrast, " ===")

  sig_dir <- file.path(OUTDIR, contrast)
  loocv_df <- read.csv(file.path(sig_dir, "selected_LOOCV_scores.csv"))
  best_cand <- read.csv(file.path(sig_dir, "selected_best_candidate.csv"))
  K_best    <- best_cand$Selected_K[1]
  alpha_best <- best_cand$Selected_alpha[1]
  observed_auc <- best_cand$AUC_LOOCV[1]

  y_true <- as.integer(loocv_df$Group == "NonSurvivor")
  preds  <- loocv_df$LOOCV_score
  n_samples <- length(y_true)

  ## - DeLong CI --
  roc_obj <- suppressWarnings(
    roc(y_true, preds, levels = c(0, 1), direction = "<", quiet = TRUE)
  )
  ci_delong <- tryCatch(
    as.numeric(ci.auc(roc_obj, method = "delong", conf.level = 0.95)),
    error = function(e) c(NA, observed_auc, NA)
  )
  ci_lo <- ci_delong[1]; ci_hi <- ci_delong[3]
  message(sprintf("  DeLong CI: AUC=%.4f [%.4f, %.4f]", observed_auc, ci_lo, ci_hi))

  ## - Bootstrap CI (as backup) --
  ci_boot <- tryCatch(
    as.numeric(ci.auc(roc_obj, method = "bootstrap", boot.n = 1000, conf.level = 0.95)),
    error = function(e) c(NA, observed_auc, NA)
  )
  message(sprintf("  Bootstrap CI: [%.4f, %.4f]", ci_boot[1], ci_boot[3]))

  ## - Permutation test ----------------------------------------
  message(sprintf("  Running %d-permutation test (fixed K=%d, alpha=%.2f)...", N_PERM, K_best, alpha_best))

  meta_sub <- as.data.frame(contrast_meta[[contrast]])
  common_ids <- intersect(colnames(cts_sym), meta_sub$Core_Sample_Id)
  Xraw <- cts_sym[, common_ids, drop = FALSE]
  meta  <- meta_sub[match(common_ids, meta_sub$Core_Sample_Id), , drop = FALSE]
  norm  <- normalize_contrast(Xraw, meta, do_batch_remove = TRUE)
  mat   <- norm$mat
  mat_z <- norm$mat_z
  y_cv  <- as.integer(meta$GroupBin == "NonSurvivor")
  group_bin <- meta$GroupBin

  null_aucs <- numeric(N_PERM)
  set.seed(42)
  for (perm_i in seq_len(N_PERM)) {
    y_perm <- sample(y_cv)
    gb_perm <- ifelse(y_perm == 1, "NonSurvivor", "Survivor")
    preds_perm <- tryCatch(
      loocv_fixed_k_alpha(mat, mat_z, y_perm, gb_perm, K_best, alpha_best),
      error = function(e) rep(0.5, length(y_perm))
    )
    null_aucs[perm_i] <- auc_roc(y_perm, preds_perm)
    if (perm_i %% 50 == 0) message(sprintf("    perm %d/%d done", perm_i, N_PERM))
  }

  perm_p <- mean(null_aucs >= observed_auc)
  message(sprintf("  Permutation p-value: %.4f (obs AUC=%.4f, null median=%.4f)",
    perm_p, observed_auc, median(null_aucs)))

  ## - Gene stability ------------------------------------------
  n <- length(y_cv)
  sig_genes_df <- read.csv(file.path(sig_dir, "selected_signature_genes.csv"))
  fixed_genes <- sig_genes_df$gene

  gene_count <- setNames(integer(length(fixed_genes)), fixed_genes)
  for (i in seq_len(n)) {
    tr <- setdiff(seq_len(n), i)
    t_tr <- limma_tstat(mat[, tr, drop = FALSE], group_bin[tr], "NonSurvivor")
    top_genes <- names(t_tr)[order(abs(t_tr), decreasing = TRUE)[seq_len(min(K_best, length(t_tr)))]]
    Mtr <- t(mat_z[top_genes, tr, drop = FALSE])
    fid <- make_inner_foldid(y_cv[tr], nfolds = min(length(tr), 5L))
    cv_fit <- suppressWarnings(cv.glmnet(
      Mtr, y_cv[tr], alpha = alpha_best, family = "binomial",
      standardize = FALSE, foldid = fid
    ))
    fit <- glmnet(Mtr, y_cv[tr], alpha = alpha_best, lambda = cv_fit$lambda.min,
                  family = "binomial", standardize = FALSE)
    coefs_i <- as.matrix(coef(fit))[-1, , drop = TRUE]
    fold_genes <- names(coefs_i)[coefs_i != 0]
    for (g in fold_genes) {
      if (g %in% names(gene_count)) gene_count[g] <- gene_count[g] + 1L
    }
  }
  stability_df <- data.frame(
    gene = names(gene_count),
    folds_selected = gene_count,
    frac_selected = gene_count / n,
    row.names = NULL
  )
  stability_df <- stability_df[order(-stability_df$frac_selected), ]

  n_stable <- sum(stability_df$frac_selected >= 0.70)
  message(sprintf("  Gene stability: %d/%d genes selected in >=70%% of folds",
    n_stable, nrow(stability_df)))

  ## - Save results -------------------------------------------
  out_dir_val <- file.path(sig_dir, "validation_ci_perm")
  dir.create(out_dir_val, recursive = TRUE, showWarnings = FALSE)

  null_df <- data.frame(perm = seq_len(N_PERM), null_auc = null_aucs)
  write.csv(null_df, file.path(out_dir_val, "permutation_null_aucs.csv"), row.names = FALSE)
  write.csv(stability_df, file.path(out_dir_val, "gene_stability.csv"), row.names = FALSE)

  png(file.path(out_dir_val, "permutation_null_histogram.png"), width = 700, height = 500, res = 96)
  hist(null_aucs, breaks = 30, col = "steelblue", border = "white",
    main = sprintf("%s\nNull AUC distribution (N=%d permutations)", contrast, N_PERM),
    xlab = "Permuted AUC", xlim = c(0, 1))
  abline(v = observed_auc, col = "red", lwd = 2.5)
  legend("topleft", legend = sprintf("Observed AUC = %.3f\np = %.4f", observed_auc, perm_p),
    col = "red", lwd = 2, bty = "n")
  dev.off()

  res_row <- data.frame(
    Contrast        = contrast,
    N               = n_samples,
    LOOCV_AUC       = observed_auc,
    DeLong_CI_lo    = ci_lo,
    DeLong_CI_hi    = ci_hi,
    Boot_CI_lo      = ci_boot[1],
    Boot_CI_hi      = ci_boot[3],
    Perm_p          = perm_p,
    N_perms         = N_PERM,
    N_stable_genes  = n_stable,
    Total_sig_genes = nrow(stability_df),
    CI_lo_gt_0.5    = !is.na(ci_lo) && ci_lo > 0.50,
    CI_lo_gt_0.6    = !is.na(ci_lo) && ci_lo > 0.60,
    Perm_p_lt_0.05  = perm_p < 0.05,
    stringsAsFactors = FALSE
  )
  results_all[[contrast]] <- res_row
  message(sprintf("  Done. CI_lo>0.6=%s  perm_p<0.05=%s  stable_genes=%d",
    res_row$CI_lo_gt_0.6, res_row$Perm_p_lt_0.05, n_stable))
}

## ---- Write summary CSV -------------------------------------
summary_df <- do.call(rbind, results_all)
write.csv(summary_df,
  file.path(OUTDIR, "validation_ci_perm_summary.csv"),
  row.names = FALSE)

message("\n=== SUMMARY ===")
print(summary_df[, c("Contrast", "N", "LOOCV_AUC", "DeLong_CI_lo", "DeLong_CI_hi",
                      "Perm_p", "CI_lo_gt_0.6", "Perm_p_lt_0.05", "N_stable_genes")])

message("Summary CSV: ", file.path(OUTDIR, "validation_ci_perm_summary.csv"))
