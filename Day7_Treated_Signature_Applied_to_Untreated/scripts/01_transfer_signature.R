## Apply the day-7 treated elastic-net signature (Fig 3C) to the untreated day-7 mice.
## Coefficients and intercept are taken from the treated fit; nothing is refit.
## Two standardizations: A) z-scored with the treated cohort's mean/SD,
## B) z-scored within the untreated cohort.

suppressPackageStartupMessages({
  .libPaths(c("./NCT_paper_pipeline/renv_lib",
              .libPaths()))
  library(readr); library(dplyr); library(matrixStats); library(edgeR); library(pROC)
})
set.seed(123)

PROJ <- "."
PIPE <- file.path(PROJ, "NCT_paper_pipeline")
OUT  <- file.path(PROJ, "Day7_Treated_Signature_Applied_to_Untreated")
RES  <- file.path(OUT, "results"); dir.create(RES, recursive = TRUE, showWarnings = FALSE)

INPUT_COUNTS <- file.path(PIPE, "01_inputs", "merged_counts.tsv")
INPUT_META   <- file.path(PIPE, "01_inputs", "merged_metadata_all.tsv")
SIG_CSV <- file.path(PIPE, "02_cache", "elastic_net_day7_day18_source_match",
                     "results_EN_signatures", "Day7_Treated_SvNS",
                     "selected_signature_genes.csv")

## ---- gene-symbol matrix (same as 10_figure3_housestyle.R) ----
strip_ver <- function(x) sub("\\.\\d+$", "", x)
map_sym <- function(e){ e <- strip_ver(e); odb <- getExportedValue("org.Mm.eg.db","org.Mm.eg.db")
  tab <- AnnotationDbi::select(odb, keys=e, keytype="ENSEMBL", columns="SYMBOL")
  o <- rep(NA_character_, length(e)); o[match(tab$ENSEMBL, e)] <- tab$SYMBOL; o }

meta_all <- read_tsv(INPUT_META, show_col_types = FALSE)
cts_tbl  <- read_tsv(INPUT_COUNTS, show_col_types = FALSE)
gene_col <- names(cts_tbl)[1]; sample_cols <- setdiff(names(cts_tbl), gene_col)
ens <- strip_ver(as.character(cts_tbl[[gene_col]]))
mat_raw <- as.matrix(cts_tbl[, sample_cols]); storage.mode(mat_raw) <- "numeric"; rownames(mat_raw) <- ens
mat_raw <- rowsum(mat_raw, rownames(mat_raw), reorder = FALSE)
sym <- map_sym(rownames(mat_raw))
keep_sym <- !is.na(sym) & nzchar(sym) & !grepl("Rik$", sym) & !grepl("^(Gm|Ps)", sym)
mat_raw <- mat_raw[keep_sym, ]; mat_sym <- rowsum(mat_raw, sym[keep_sym], reorder = FALSE)

## ---- Day 7 metadata: treated (training cohort) + untreated (target cohort) ----
m_day7 <- meta_all %>%
  filter(Timepoint == "Day7", Survival_Status %in% c("survivor", "non-survivor")) %>%
  mutate(Treatment = tolower(Treatment),
         GroupBin  = ifelse(tolower(Survival_Status) == "non-survivor", "NonSurvivor", "Survivor")) %>%
  as.data.frame()
common <- intersect(colnames(mat_sym), m_day7$Core_Sample_Id)
m_day7 <- m_day7 %>% filter(Core_Sample_Id %in% common) %>%
  arrange(match(Core_Sample_Id, common)) %>% as.data.frame()
stopifnot(length(unique(m_day7$Batch)) == 1)   # all Day7 is batch_2 -> no batch removal

is_tr <- m_day7$Treatment == "y"
message("Day 7 samples: ", nrow(m_day7), "  treated=", sum(is_tr), "  untreated=", sum(!is_tr))
message("  treated   ", sum(is_tr & m_day7$GroupBin=="Survivor"), "S / ",
        sum(is_tr & m_day7$GroupBin=="NonSurvivor"), "NS")
message("  untreated ", sum(!is_tr & m_day7$GroupBin=="Survivor"), "S / ",
        sum(!is_tr & m_day7$GroupBin=="NonSurvivor"), "NS")

## ---- ONE joint normalization over all 23 Day 7 samples ----
X <- mat_sym[, m_day7$Core_Sample_Id, drop = FALSE]
d <- DGEList(X); d <- calcNormFactors(d, "TMM")
keep <- rowSums(cpm(d, log = FALSE) >= 1) >= ceiling(0.25 * ncol(d))
d <- d[keep, , keep.lib.sizes = FALSE]
logcpm <- cpm(d, log = TRUE, prior.count = 1)
logcpm <- logcpm[rowVars(logcpm) > 0, ]
message("Genes after joint Day7 filter: ", nrow(logcpm))

## ---- frozen signature ----
sig <- read_csv(SIG_CSV, show_col_types = FALSE) %>% arrange(desc(abs(coef)))
missing <- setdiff(sig$gene, rownames(logcpm))
if (length(missing)) stop("Signature genes lost in the joint Day7 filter: ", paste(missing, collapse=", "))
genes <- sig$gene; b <- sig$coef
message("Signature genes present: ", length(genes), "/", nrow(sig))

## ---- two standardizations, both applied to the frozen coefficients ----
# z-score every signature gene using the mean/SD of a chosen reference set of
# columns, then return the linear predictor for all 23 samples.
eta_using_reference <- function(ref_cols) {
  M  <- logcpm[genes, , drop = FALSE]
  mu <- rowMeans(M[, ref_cols, drop = FALSE])
  sd <- rowSds(M[, ref_cols, drop = FALSE]); sd[sd == 0] <- 1
  Z  <- (M - mu) / sd
  as.numeric(t(Z) %*% b)
}
eta_A <- eta_using_reference(is_tr)    # strict transfer  (treated reference)
eta_B <- eta_using_reference(!is_tr)   # re-standardized  (untreated reference)

## Intercept fit once on the treated cohort with the linear predictor as offset
## (same construction as 10_figure3_housestyle.R:78); reused for both variants.
y_tr <- as.integer(m_day7$GroupBin[is_tr] == "NonSurvivor")
b0 <- suppressWarnings(coef(glm(y_tr ~ 1 + offset(eta_A[is_tr]), family = binomial))[1])
message(sprintf("Treated-cohort intercept b0 = %.4f", b0))

## ---- check: does the joint normalization change the model's footing? ----
## Rebuild the treated cohort exactly as the published figure does (TMM + CPM
## filter over the 12 treated samples only, then z-score) and compare its linear
## predictor to the joint-space one. A near-1 correlation means the joint
## normalization is innocuous and any transfer failure is not an artifact of it.
auc_of <- function(y, s) suppressWarnings(as.numeric(auc(roc(y, s, levels = c(0,1),
                                                             direction = "<", quiet = TRUE))))
auc_treated_joint <- auc_of(y_tr, eta_A[is_tr])

d_tr <- DGEList(mat_sym[, m_day7$Core_Sample_Id[is_tr], drop = FALSE])
d_tr <- calcNormFactors(d_tr, "TMM")
keep_tr <- rowSums(cpm(d_tr, log = FALSE) >= 1) >= ceiling(0.25 * ncol(d_tr))
d_tr <- d_tr[keep_tr, , keep.lib.sizes = FALSE]
M_tr <- cpm(d_tr, log = TRUE, prior.count = 1)[genes, , drop = FALSE]
Z_tr <- (M_tr - rowMeans(M_tr)) / ifelse(rowSds(M_tr) == 0, 1, rowSds(M_tr))
eta_treated_only <- as.numeric(t(Z_tr) %*% b)
cor_joint_vs_own <- cor(eta_treated_only, eta_A[is_tr])
auc_treated_own  <- auc_of(y_tr, eta_treated_only)
message(sprintf("Sanity: treated eta joint vs treated-only normalization r = %.4f (apparent AUC %.3f vs %.3f)",
                cor_joint_vs_own, auc_treated_joint, auc_treated_own))

## ---- untreated cohort: scores vs actual outcome ----
mu <- m_day7[!is_tr, ]
y_un <- as.integer(mu$GroupBin == "NonSurvivor")
sc <- data.frame(
  Core_Sample_Id = mu$Core_Sample_Id, Sample_Name = mu$Sample_Name,
  Actual_Group = factor(mu$GroupBin, levels = c("Survivor", "NonSurvivor")),
  eta_A = eta_A[!is_tr], prob_A = plogis(b0 + eta_A[!is_tr]),
  eta_B = eta_B[!is_tr], prob_B = plogis(b0 + eta_B[!is_tr]),
  stringsAsFactors = FALSE)

perf <- function(eta, tag) {
  r  <- suppressWarnings(roc(y_un, eta, levels = c(0,1), direction = "<", quiet = TRUE))
  ci <- suppressWarnings(ci.auc(r, method = "delong"))
  w  <- suppressWarnings(wilcox.test(eta ~ y_un))
  data.frame(variant = tag, n = length(y_un),
             AUC_vs_actual = as.numeric(auc(r)),
             CI_low = as.numeric(ci)[1], CI_high = as.numeric(ci)[3],
             wilcox_p = unname(w$p.value),
             mean_survivor = mean(eta[y_un == 0]), mean_nonsurvivor = mean(eta[y_un == 1]))
}
perf_tbl <- rbind(perf(sc$eta_A, "A_strict_transfer"), perf(sc$eta_B, "B_restandardized"))
print(perf_tbl)

## ---- predicted groups (median split of the strict-transfer score) ----
## Median split, not a 0.5 probability cut: the intercept is calibrated to a
## different cohort, so a fixed probability threshold can easily yield a 1-vs-10
## split that DESeq2 cannot fit. Median split guarantees 5/6 at n=11.
med <- median(sc$eta_A)
sc$Predicted_Group_medianA <- ifelse(sc$eta_A >= med, "Predicted_NonSurvivor", "Predicted_Survivor")
sc$Predicted_Group_p50A    <- ifelse(sc$prob_A > 0.5, "Predicted_NonSurvivor", "Predicted_Survivor")

## Margin exclusion: keep the median cut but drop the N_MARGIN_DROP mice sitting
## closest to it, whose side is essentially decided by noise (the two nearest are
## 0.04 apart). Widens the gap between groups at the cost of n and of removing
## the hardest-to-classify samples - which inflates apparent separation even
## under the null, so this contrast is descriptive, not a test.
N_MARGIN_DROP <- 2
margin_rank <- rank(abs(sc$eta_A - med), ties.method = "first")
sc$Predicted_Group_margin <- ifelse(margin_rank <= N_MARGIN_DROP, NA_character_,
                                    sc$Predicted_Group_medianA)

## Continuous: the score standardized to SD = 1, for use as a numeric covariate
## in DESeq2. No threshold at all - uses every mouse and the full score.
sc$eta_A_z <- as.numeric(scale(sc$eta_A))

sc <- sc[order(sc$eta_A), ]

message("Predicted (median split): ",
        sum(sc$Predicted_Group_medianA == "Predicted_Survivor"), " pred-Survivor / ",
        sum(sc$Predicted_Group_medianA == "Predicted_NonSurvivor"), " pred-NonSurvivor")
message("Predicted (p>0.5):        ",
        sum(sc$Predicted_Group_p50A == "Predicted_Survivor"), " pred-Survivor / ",
        sum(sc$Predicted_Group_p50A == "Predicted_NonSurvivor"), " pred-NonSurvivor")
message("Predicted (margin-excl):  ",
        sum(sc$Predicted_Group_margin == "Predicted_Survivor", na.rm = TRUE), " pred-Survivor / ",
        sum(sc$Predicted_Group_margin == "Predicted_NonSurvivor", na.rm = TRUE), " pred-NonSurvivor",
        "  (dropped: ", paste(sc$Sample_Name[is.na(sc$Predicted_Group_margin)], collapse = ", "), ")")

conf <- bind_rows(
  as.data.frame(table(Predicted = sc$Predicted_Group_medianA, Actual = sc$Actual_Group)) %>%
    mutate(rule = "median_split_etaA"),
  as.data.frame(table(Predicted = sc$Predicted_Group_p50A, Actual = sc$Actual_Group)) %>%
    mutate(rule = "prob_gt_0.5"))
concord_med <- mean((sc$Predicted_Group_medianA == "Predicted_NonSurvivor") == (sc$Actual_Group == "NonSurvivor"))
concord_p50 <- mean((sc$Predicted_Group_p50A    == "Predicted_NonSurvivor") == (sc$Actual_Group == "NonSurvivor"))
message(sprintf("Concordance with actual outcome: median split %.3f, p>0.5 %.3f", concord_med, concord_p50))

## ---- treated-cohort scores kept for the diagnostics figure ----
treated_sc <- data.frame(
  Core_Sample_Id = m_day7$Core_Sample_Id[is_tr], Sample_Name = m_day7$Sample_Name[is_tr],
  Cohort = "Treated (training)",
  Actual_Group = factor(m_day7$GroupBin[is_tr], levels = c("Survivor", "NonSurvivor")),
  eta_A = eta_A[is_tr], prob_A = plogis(b0 + eta_A[is_tr]))

summary_tbl <- data.frame(
  item = c("n_treated_training", "n_untreated_target",
           "n_untreated_survivor", "n_untreated_nonsurvivor",
           "n_signature_genes", "genes_after_joint_filter",
           "intercept_b0_from_treated",
           "AUC_treated_apparent_joint_space", "AUC_treated_apparent_own_space",
           "cor_eta_joint_vs_treated_only_norm", "AUC_published_treated_LOOCV",
           "AUC_transfer_A_vs_actual", "AUC_transfer_B_vs_actual",
           "wilcox_p_A", "wilcox_p_B",
           "concordance_medianA_vs_actual", "concordance_p50A_vs_actual",
           "n_pred_survivor_medianA", "n_pred_nonsurvivor_medianA",
           "n_pred_survivor_margin", "n_pred_nonsurvivor_margin", "n_dropped_margin"),
  value = c(sum(is_tr), sum(!is_tr),
            sum(y_un == 0), sum(y_un == 1),
            length(genes), nrow(logcpm),
            round(b0, 4),
            round(auc_treated_joint, 4), round(auc_treated_own, 4),
            round(cor_joint_vs_own, 4), 0.8056,
            round(perf_tbl$AUC_vs_actual[1], 4), round(perf_tbl$AUC_vs_actual[2], 4),
            signif(perf_tbl$wilcox_p[1], 4), signif(perf_tbl$wilcox_p[2], 4),
            round(concord_med, 4), round(concord_p50, 4),
            sum(sc$Predicted_Group_medianA == "Predicted_Survivor"),
            sum(sc$Predicted_Group_medianA == "Predicted_NonSurvivor"),
            sum(sc$Predicted_Group_margin == "Predicted_Survivor", na.rm = TRUE),
            sum(sc$Predicted_Group_margin == "Predicted_NonSurvivor", na.rm = TRUE),
            N_MARGIN_DROP),
  stringsAsFactors = FALSE)

write_csv(sc,          file.path(RES, "untreated_day7_en_scores.csv"))
write_csv(treated_sc,  file.path(RES, "treated_day7_en_scores.csv"))
write_csv(conf,        file.path(RES, "predicted_vs_actual_confusion.csv"))
write_csv(perf_tbl,    file.path(RES, "transfer_performance.csv"))
write_csv(summary_tbl, file.path(RES, "contrast_summary.csv"))

message("\n=== untreated Day 7 scores ===")
print(sc[, c("Sample_Name", "Actual_Group", "eta_A", "prob_A", "eta_B",
             "Predicted_Group_medianA")], row.names = FALSE)
message("\nTreated apparent AUC in joint space: ", round(auc_treated_joint, 4),
        "  (published LOOCV 0.806)")
message("Written to ", RES)
