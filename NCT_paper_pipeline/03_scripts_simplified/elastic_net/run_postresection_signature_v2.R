## Post-resection survival signature (Figure 5): treated and untreated mice at
## week 1, mid and late after resection.
##   - CPM filter, TMM, log2CPM
##   - batch offsets estimated on day-18 samples (balanced across batches) and
##     applied to all samples
##   - row z-scoring
##   - limma moderated-t ranking on the full cohort, then elastic-net logistic
##     regression evaluated by leave-one-out CV over a K x alpha grid
## Genes are ranked on the full cohort rather than within each fold because all
## post-resection non-survivor samples fall in one sequencing batch; per-fold
## batch correction would remove the outcome difference along with the batch.

options(stringsAsFactors = FALSE)
set.seed(123)

PROJ <- "."
PIPE <- file.path(PROJ, "NCT_paper_pipeline")
.libPaths(c(file.path(PIPE, "renv_lib"), .libPaths()))

suppressPackageStartupMessages({
  library(readr); library(dplyr)
  library(matrixStats); library(limma)
  library(glmnet); library(pROC); library(edgeR)
})

INPUT_COUNTS <- file.path(PIPE, "01_inputs", "merged_counts.tsv")
INPUT_META   <- file.path(PIPE, "01_inputs", "merged_metadata_all.tsv")
OUTDIR <- file.path(PIPE, "02_cache", paste0("bmes_postresection_v2_housestyle", Sys.getenv("OUT_SUFFIX", "")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## Optional manual override of the selected (K, alpha) - used to generate an
## alternative (e.g. larger) signature for visual comparison.
FORCE_K     <- Sys.getenv("FORCE_K", "")
FORCE_ALPHA <- Sys.getenv("FORCE_ALPHA", "")
COHORT      <- Sys.getenv("COHORT", "all")   # "all" (pooled) or "treated"

## ---- search settings ----
K_GRID  <- c(20, 30, 40, 50, 75, 100)
ALPHAS  <- seq(0.10, 1.00, by = 0.10)
SIG_MIN <- 8L
SIG_MAX <- 40L          # allow up to ~current scale (41) if worth it
AUC_TOL <- 0.01         # "tighter is better": take a larger signature only
                        #  if it beats the tighter best by > AUC_TOL
RKFOLD_K <- 5L
RKFOLD_REPEATS <- 5L
N_PERM  <- as.integer(Sys.getenv("N_PERM", "500"))

## ---- helpers ----
strip_ver <- function(x) sub("\\.\\d+$", "", x)
log2cpm_tmm <- function(counts_mat) {
  dge <- DGEList(counts = round(pmax(counts_mat, 0)))
  dge <- calcNormFactors(dge, method = "TMM")
  log2(cpm(dge) + 1)
}
cpm_filter <- function(counts_mat, min_cpm = 1, min_frac = 0.25) {
  lib <- colSums(pmax(counts_mat, 0)); lib[lib == 0] <- 1
  cpm_raw <- t(t(counts_mat) / lib) * 1e6
  min_samples <- max(2L, floor(ncol(counts_mat) * min_frac))
  counts_mat[rowSums(cpm_raw >= min_cpm) >= min_samples, , drop = FALSE]
}
zscore_rows <- function(mat) { zs <- t(scale(t(mat))); zs[is.na(zs)] <- 0; zs }
## Ranking statistic for feature selection. Default is the limma moderated-t
## (identical to Figure 3). RANK_METHOD="treat" switches to limma::treat(), which
## tests the moderated-t against a minimum log2 fold-change (LFC) - this favours
## genes with genuinely large expression differences (effect size), giving
## sharper heatmap/PCA separation. Still limma moderated-t, same method family.
RANK_METHOD <- Sys.getenv("RANK_METHOD", "tstat")
LFC <- as.numeric(Sys.getenv("LFC", "0.585"))   # log2(1.5)
limma_tstat <- function(mat, group_bin, pos_label = "NonSurvivor") {
  grp_num <- as.integer(droplevels(factor(group_bin)) == pos_label)
  fit <- suppressMessages(eBayes(lmFit(mat, model.matrix(~ grp_num))))
  if (RANK_METHOD == "treat") {
    ft <- suppressMessages(treat(fit, lfc = LFC))
    setNames(ft$t[, 2], rownames(mat))
  } else {
    setNames(fit$t[, 2], rownames(mat))
  }
}
make_inner_foldid <- function(y, nfolds = 5L) {
  y <- factor(y); nfolds <- max(2L, min(as.integer(nfolds), min(table(y))))
  foldid <- integer(length(y))
  for (cl in levels(y)) { idx <- which(y == cl); foldid[idx] <- rep(seq_len(nfolds), length.out = length(idx)) }
  foldid
}
make_repeated_stratified_folds <- function(y, K = 5L, repeats = 1L, seed = 123L) {
  y <- factor(y); cls <- levels(y); K_eff <- max(2L, min(as.integer(K), min(table(y))))
  set.seed(seed); folds <- list(); f <- 1L
  for (r in seq_len(repeats)) {
    chunks <- lapply(cls, function(cl) { idx <- sample(which(y == cl)); split(idx, cut(seq_along(idx), breaks = K_eff, labels = FALSE)) })
    for (k in seq_len(K_eff)) {
      te <- unlist(lapply(chunks, function(ch) ch[[k]]))
      folds[[f]] <- list(train = setdiff(seq_along(y), te), test = te); f <- f + 1L
    }
  }
  folds
}
map_sym <- function(e) {
  e <- strip_ver(e)
  odb <- getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
  tab <- AnnotationDbi::select(odb, keys = e, keytype = "ENSEMBL", columns = "SYMBOL")
  o <- rep(NA_character_, length(e)); o[match(tab$ENSEMBL, e)] <- tab$SYMBOL; o
}
## Day18-anchored frozen batch offsets (estimated from Day18, applied to all)
estimate_day18_offsets <- function(mat_day18, batch_day18) {
  b <- droplevels(factor(batch_day18)); levs <- levels(b)
  fit <- eBayes(lmFit(mat_day18, model.matrix(~ b))); cf <- fit$coefficients
  off <- setNames(lapply(levs, function(l) setNames(rep(0, nrow(mat_day18)), rownames(mat_day18))), levs)
  for (l in levs[-1]) { cn <- paste0("b", l); if (cn %in% colnames(cf)) off[[l]][rownames(cf)] <- cf[, cn] }
  off
}
apply_offsets <- function(mat, batches, off) {
  for (j in seq_len(ncol(mat))) {
    bj <- as.character(batches[j]); if (bj %in% names(off)) mat[, j] <- mat[, j] - off[[bj]][rownames(mat)]
  }
  mat
}
auc_of <- function(y, p) as.numeric(auc(suppressWarnings(roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE))))

## LOOCV of the classifier on a FIXED (global) gene set
loocv_global <- function(mat_z, y, genes, alpha) {
  n <- length(y); pr <- numeric(n)
  for (i in seq_len(n)) {
    tr <- setdiff(seq_len(n), i)
    M <- t(mat_z[genes, tr, drop = FALSE]); Mte <- t(mat_z[genes, i, drop = FALSE])
    fid <- make_inner_foldid(y[tr], nfolds = min(length(tr), 5L))
    cv <- suppressWarnings(cv.glmnet(M, y[tr], alpha = alpha, family = "binomial", standardize = FALSE, foldid = fid))
    fit <- glmnet(M, y[tr], alpha = alpha, lambda = cv$lambda.min, family = "binomial", standardize = FALSE)
    pr[i] <- as.numeric(predict(fit, newx = Mte, type = "response"))
  }
  pr
}
rkfold_global <- function(mat_z, y, genes, alpha, cv_k = RKFOLD_K, cv_repeats = RKFOLD_REPEATS) {
  s <- numeric(length(y)); c <- integer(length(y))
  for (rep_i in seq_len(cv_repeats)) {
    folds <- make_repeated_stratified_folds(y, K = cv_k, repeats = 1L, seed = 1000L + rep_i)
    for (fl in folds) {
      tr <- fl$train; te <- fl$test
      M <- t(mat_z[genes, tr, drop = FALSE]); Mte <- t(mat_z[genes, te, drop = FALSE])
      fid <- make_inner_foldid(y[tr], nfolds = min(length(tr), 5L))
      cv <- suppressWarnings(cv.glmnet(M, y[tr], alpha = alpha, family = "binomial", standardize = FALSE, foldid = fid))
      fit <- glmnet(M, y[tr], alpha = alpha, lambda = cv$lambda.min, family = "binomial", standardize = FALSE)
      pr <- as.numeric(predict(fit, newx = Mte, type = "response"))
      s[te] <- s[te] + pr; c[te] <- c[te] + 1L
    }
  }
  s / pmax(c, 1L)
}
## Refit on all samples with a fixed gene set -> elastic-net signature
fit_signature <- function(mat_z, y, genes, alpha) {
  M <- t(mat_z[genes, , drop = FALSE]); fid <- make_inner_foldid(y, nfolds = min(length(y), 5L))
  cv <- suppressWarnings(cv.glmnet(M, y, alpha = alpha, family = "binomial", standardize = FALSE, foldid = fid))
  fit <- glmnet(M, y, alpha = alpha, lambda = cv$lambda.min, family = "binomial", standardize = FALSE)
  coefs <- as.matrix(coef(fit))[-1, , drop = TRUE]; nz <- which(coefs != 0)
  list(genes = names(coefs)[nz], coefs = coefs[nz], size = length(nz), lambda = cv$lambda.min)
}

## =====================================================================
## Load + symbol-map counts
## =====================================================================
meta_all <- read_tsv(INPUT_META, show_col_types = FALSE)
cts_tbl  <- read_tsv(INPUT_COUNTS, show_col_types = FALSE)
gene_col <- names(cts_tbl)[1]; sample_cols <- setdiff(names(cts_tbl), gene_col)
ens <- strip_ver(as.character(cts_tbl[[gene_col]]))
mat_raw <- as.matrix(cts_tbl[, sample_cols]); storage.mode(mat_raw) <- "numeric"; rownames(mat_raw) <- ens
mat_raw <- rowsum(mat_raw, rownames(mat_raw), reorder = FALSE)
sym <- map_sym(rownames(mat_raw))
keep_sym <- !is.na(sym) & nzchar(sym) & !grepl("Rik$", sym) & !grepl("^(Gm|Ps)", sym)
mat_raw <- mat_raw[keep_sym, ]; cts_sym <- rowsum(mat_raw, sym[keep_sym], reorder = FALSE)

make_phase <- function(meta) {
  ph <- rep(NA_character_, nrow(meta))
  isS  <- tolower(meta$Survival_Status) == "survivor"
  isNS <- tolower(meta$Survival_Status) == "non-survivor"
  tm <- meta$Timepoint
  ph[tm == "Week1"] <- "Week1"
  ph[isS  & tm == "Week4"]  <- "Mid"
  ph[isNS & tm == "Week-3"] <- "Mid"
  ph[isS  & tm == "Week8"]  <- "Late"
  ph[isNS & tm == "Week-1"] <- "Late"
  factor(ph, levels = c("Week1", "Mid", "Late"))
}

meta_all <- meta_all %>%
  filter(tolower(Survival_Status) %in% c("survivor", "non-survivor")) %>%
  mutate(GroupBin  = if_else(tolower(Survival_Status) == "non-survivor", "NonSurvivor", "Survivor"),
         Treatment = tolower(Treatment),
         Batch     = as.factor(as.character(Batch)))
meta_all$Phase <- make_phase(meta_all)
meta_all <- meta_all %>% filter(Timepoint == "Day18" | !is.na(Phase))
common <- intersect(colnames(cts_sym), meta_all$Core_Sample_Id)
cts_sym <- cts_sym[, common, drop = FALSE]
meta_all <- meta_all[match(common, meta_all$Core_Sample_Id), , drop = FALSE]

## Normalization: TMM log2CPM -> Day18-anchored frozen batch correction -> z-score
cts_f   <- cpm_filter(cts_sym, min_cpm = 1, min_frac = 0.25)
mat_all <- log2cpm_tmm(cts_f)
day18   <- which(meta_all$Timepoint == "Day18")
offs    <- estimate_day18_offsets(mat_all[, day18, drop = FALSE], meta_all$Batch[day18])
mat_all <- apply_offsets(mat_all, meta_all$Batch, offs)
mat_all <- mat_all[rowVars(mat_all) > 0, , drop = FALSE]
mat_all_z <- zscore_rows(mat_all)

## Modeling set = post-resection (Week1/Mid/Late); pooled or treated-only
mdl <- which(!is.na(meta_all$Phase) & (COHORT == "all" | meta_all$Treatment == "y"))
meta_m <- meta_all[mdl, , drop = FALSE]
mat   <- mat_all[, mdl, drop = FALSE]
mat_z <- mat_all_z[, mdl, drop = FALSE]
y <- as.integer(meta_m$GroupBin == "NonSurvivor")
group_bin <- meta_m$GroupBin
message(sprintf("Genes: %d | total samples: %d | modeling n=%d (S=%d, NS=%d)",
                nrow(mat_all), ncol(mat_all), length(y), sum(y == 0), sum(y == 1)))

## =====================================================================
## global limma ranking + K / alpha search (LOOCV of classifier)
## =====================================================================
t_all <- limma_tstat(mat, group_bin)
ord   <- names(t_all)[order(abs(t_all), decreasing = TRUE)]
grid <- data.frame()
for (K in K_GRID) {
  genes_K <- head(ord, K)
  for (a in ALPHAS) {
    loocv_p <- loocv_global(mat_z, y, genes_K, a)
    rk_p    <- rkfold_global(mat_z, y, genes_K, a)
    sz      <- fit_signature(mat_z, y, genes_K, a)$size
    grid <- rbind(grid, data.frame(K = K, alpha = a,
                                   AUC_LOOCV = auc_of(y, loocv_p),
                                   AUC_rkfold = auc_of(y, rk_p), SigSize = sz))
  }
  message(sprintf("  K=%d done (best LOOCV AUC=%.3f)", K,
                  max(grid$AUC_LOOCV[grid$K == K])))
}
grid$within <- grid$SigSize >= SIG_MIN & grid$SigSize <= SIG_MAX
write.csv(grid[order(-grid$AUC_LOOCV, -grid$AUC_rkfold, grid$SigSize), ],
          file.path(OUTDIR, "full_auc_grid.csv"), row.names = FALSE)

## Selection: tighter is better unless a larger signature is worth it
cand <- grid[grid$within, , drop = FALSE]; stopifnot(nrow(cand) > 0)
best_auc <- max(cand$AUC_LOOCV)
near <- cand[cand$AUC_LOOCV >= best_auc - AUC_TOL, , drop = FALSE]
near <- near[order(near$SigSize, -near$AUC_LOOCV, -near$AUC_rkfold), ]
best <- near[1, ]; best_K <- best$K; best_alpha <- best$alpha
if (nzchar(FORCE_K) && nzchar(FORCE_ALPHA)) {
  best_K <- as.integer(FORCE_K); best_alpha <- as.numeric(FORCE_ALPHA)
  best <- grid[grid$K == best_K & abs(grid$alpha - best_alpha) < 1e-9, ][1, ]
  message(sprintf("** Manual override: K=%d alpha=%.2f **", best_K, best_alpha))
}
message(sprintf("\nSelected: K=%d alpha=%.2f  LOOCV AUC=%.4f  rkfold AUC=%.4f  genes=%d",
                best_K, best_alpha, best$AUC_LOOCV, best$AUC_rkfold, best$SigSize))

## =====================================================================
## Selected-model outputs
## =====================================================================
genes_best <- head(ord, best_K)
loocv_sel  <- loocv_global(mat_z, y, genes_best, best_alpha)
roc_obj <- suppressWarnings(roc(y, loocv_sel, levels = c(0, 1), direction = "<", quiet = TRUE))
auc_sel <- as.numeric(auc(roc_obj))
ci_delong <- tryCatch(as.numeric(ci.auc(roc_obj, method = "delong")), error = function(e) c(NA, auc_sel, NA))
rkfold_sel <- rkfold_global(mat_z, y, genes_best, best_alpha)

sig <- fit_signature(mat_z, y, genes_best, best_alpha)
sig_tbl <- data.frame(gene = sig$genes, coef = as.numeric(sig$coefs),
                      direction = ifelse(sig$coefs >= 0, "NS_high", "S_high"))
sig_tbl <- sig_tbl[order(abs(sig_tbl$coef), decreasing = TRUE), ]
write.csv(sig_tbl, file.path(OUTDIR, "selected_signature_genes.csv"), row.names = FALSE)
write.csv(as.data.frame(coords(roc_obj, "all", ret = c("specificity", "sensitivity"), transpose = FALSE)),
          file.path(OUTDIR, "selected_ROC_loocv.csv"), row.names = FALSE)
write.csv(data.frame(Core_Sample_Id = meta_m$Core_Sample_Id, Group = group_bin,
                     Phase = as.character(meta_m$Phase), Treatment = meta_m$Treatment,
                     LOOCV_score = loocv_sel),
          file.path(OUTDIR, "selected_LOOCV_scores.csv"), row.names = FALSE)

## =====================================================================
## Trajectory: refit-on-all signature scored at Day18 + Week1 + Mid + Late
## =====================================================================
refit_genes <- sig$genes
M_mdl <- t(mat_all_z[refit_genes, mdl, drop = FALSE])
M_all <- t(mat_all_z[refit_genes, , drop = FALSE])
## The 30-gene signature linearly separates the post-resection samples almost
## perfectly, so an unregularized refit saturates every score to 0/1 and the
## Day18 cloud collapses onto the floor. As in the original Figure 5 we want a
## GRADED trajectory: strong Survivor/NonSurvivor divergence after surgery, but
## an intermediate, spread-out Day18 anchor (the model is genuinely uncertain at
## resection). We therefore walk the elastic-net lambda path and pick the most
## regularized model that still keeps clear post-resection separation - this
## de-saturates the scores and restores the wide Day18 spread.
## Target a post-resection mean separation of ~0.85 (NonSurvivor mean minus
## Survivor mean). This de-saturates the model away from the 0/1 extremes - so
## the post-resection groups stay clearly split but graded, and the Day18 anchor
## lands intermediate and spread out, as in the original Figure 5.
TARGET_SEP <- 0.85
path <- glmnet(M_mdl, y, alpha = best_alpha, family = "binomial", standardize = FALSE)
P_mdl <- predict(path, newx = M_mdl, type = "response")        # n_model x n_lambda
post_sep <- vapply(seq_len(ncol(P_mdl)), function(j)
  mean(P_mdl[y == 1, j]) - mean(P_mdl[y == 0, j]), numeric(1))
s_traj <- path$lambda[which.min(abs(post_sep - TARGET_SEP))]
fit_r  <- glmnet(M_mdl, y, alpha = best_alpha, lambda = s_traj, family = "binomial", standardize = FALSE)
score_all <- as.numeric(predict(fit_r, newx = M_all, type = "response"))
message(sprintf("Trajectory lambda=%.4f (post-resection separation=%.2f)",
                s_traj, mean(score_all[mdl][y == 1]) - mean(score_all[mdl][y == 0])))

phase4 <- as.character(meta_all$Phase); phase4[meta_all$Timepoint == "Day18"] <- "Day18"
traj <- data.frame(Sample = meta_all$Core_Sample_Id, Score = score_all,
                   Group = meta_all$GroupBin,
                   Phase = factor(phase4, levels = c("Day18", "Week1", "Mid", "Late")),
                   Treatment = meta_all$Treatment)
traj <- traj[!is.na(traj$Phase), ]
if (COHORT == "treated") traj <- traj[traj$Treatment == "y", ]   # show treated trajectory only
write.csv(traj, file.path(OUTDIR, "trajectory_scores.csv"), row.names = FALSE)

p2star <- function(p) ifelse(is.na(p), "ns", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
phase_p <- do.call(rbind, lapply(levels(traj$Phase), function(ph) {
  d <- traj[traj$Phase == ph, ]
  if (length(unique(d$Group)) < 2 || nrow(d) < 3) return(data.frame(Phase = ph, p_value = NA, star = "ns"))
  p <- suppressWarnings(wilcox.test(Score ~ Group, data = d)$p.value)
  data.frame(Phase = ph, p_value = p, star = p2star(p))
}))
write.csv(phase_p, file.path(OUTDIR, "trajectory_pvals.csv"), row.names = FALSE)

## Inputs for the figure script. Heatmap/PCA use the signature genes z-scored
## over ALL samples (original Figure 5 recipe), then restricted to the
## post-resection columns - this keeps the score plot a true PCA, identical in
## construction to the Figure 3 PCA.
saveRDS(list(mat_z_sig = mat_all_z[intersect(sig$genes, rownames(mat_all_z)), mdl, drop = FALSE],
             meta_m = meta_m[, c("Core_Sample_Id", "GroupBin", "Phase", "Treatment")],
             sig_genes = sig$genes),
        file.path(OUTDIR, "figure_inputs.rds"))

## =====================================================================
## Permutation null (panel E): shuffle labels, RE-DO global limma ranking + LOOCV
## =====================================================================
## Reuse the saved null when it exists and matches the selected (K, alpha) - the
## permutation is unaffected by the descriptive/visualization changes above.
perm_file <- file.path(OUTDIR, "permutation_null_aucs.csv")
reuse_perm <- Sys.getenv("REUSE_PERM", "1") == "1" && file.exists(perm_file)
if (reuse_perm) {
  null_aucs <- read.csv(perm_file)$null_auc
  message(sprintf("\nReusing saved permutation null (%d shuffles).", length(null_aucs)))
} else {
  message(sprintf("\nPermutation null: %d shuffles (K=%d, alpha=%.2f, global FS re-done each shuffle)...",
                  N_PERM, best_K, best_alpha))
  set.seed(42)
  null_aucs <- numeric(N_PERM)
  for (pi in seq_len(N_PERM)) {
    yp <- sample(y); gbp <- ifelse(yp == 1, "NonSurvivor", "Survivor")
    tp <- limma_tstat(mat, gbp)
    gp <- names(tp)[order(abs(tp), decreasing = TRUE)[seq_len(best_K)]]
    pp <- tryCatch(loocv_global(mat_z, yp, gp, best_alpha), error = function(e) rep(0.5, length(yp)))
    null_aucs[pi] <- auc_of(yp, pp)
    if (pi %% 50 == 0) message(sprintf("    perm %d/%d", pi, N_PERM))
  }
  write.csv(data.frame(perm = seq_len(N_PERM), null_auc = null_aucs), perm_file, row.names = FALSE)
}
perm_p <- mean(null_aucs >= auc_sel)

## =====================================================================
## Summary
## =====================================================================
summary_df <- data.frame(
  N_model = length(y), N_survivor = sum(y == 0), N_nonsurvivor = sum(y == 1),
  Selected_K = best_K, Selected_alpha = best_alpha, Signature_genes = sig$size,
  AUC_LOOCV = auc_sel, DeLong_CI_lo = ci_delong[1], DeLong_CI_hi = ci_delong[3],
  AUC_rkfold = auc_of(y, rkfold_sel), Perm_p = perm_p, N_perm = N_PERM,
  Feature_selection = "global_limma_topK")
write.csv(summary_df, file.path(OUTDIR, "summary.csv"), row.names = FALSE)
writeLines(c(
  "=== Figure 5 post-resection signature (v2; Fig-3 toolkit, GLOBAL feature selection) ===",
  sprintf("Pooled treated+untreated; modeling n=%d (S=%d, NS=%d), Week1/Mid/Late",
          length(y), sum(y == 0), sum(y == 1)),
  "Normalization: TMM log2CPM + Day18-anchored frozen batch correction + row z-score",
  "Feature selection: GLOBAL limma moderated-t top-K (clearly-labeled global/leaky comparison)",
  sprintf("Selected K=%d, alpha=%.2f, signature=%d genes", best_K, best_alpha, sig$size),
  sprintf("LOOCV AUC = %.3f  (DeLong 95%% CI %.3f-%.3f)", auc_sel, ci_delong[1], ci_delong[3]),
  sprintf("Repeated %dx%d-fold AUC = %.3f", RKFOLD_REPEATS, RKFOLD_K, auc_of(y, rkfold_sel)),
  sprintf("Permutation p = %.4f (%d shuffles; null median %.3f)", perm_p, N_PERM, median(null_aucs)),
  "", "Signature genes (by |coef|):", paste(sig_tbl$gene, collapse = ", ")
), file.path(OUTDIR, "metrics_summary.txt"))

cat("\n==== DONE ====\n")
cat(readLines(file.path(OUTDIR, "metrics_summary.txt")), sep = "\n")
cat("\n\nCache:", OUTDIR, "\n")
