## Compares GSEA results for the predicted split in untreated mice with the observed
## survival contrasts in untreated and treated mice (cosine similarity of NES).

options(stringsAsFactors = FALSE)
.libPaths(c("./NCT_paper_pipeline/renv_lib", .libPaths()))
suppressPackageStartupMessages({ library(readr); library(dplyr) })

PROJ <- "."
SRC  <- file.path(PROJ, "Day7_Treated_Signature_Applied_to_Untreated")
REB  <- file.path(PROJ, "Manuscript_Revision_2026-07", "gsea_rebuild")
OUT  <- file.path(PROJ, "Supplementary_Figure_Day7_Transfer")

source(file.path(REB, "pathway_sets.R"))   # read-only: locks the 31-pathway catalog + module order

## ---- assert the split is the one the GSEA was run on ----
calls <- read_csv(file.path(OUT, "tables", "S_per_mouse_scores_and_calls.csv"), show_col_types = FALSE) %>%
  filter(Arm == "Untreated")
old <- read_csv(file.path(SRC, "results", "untreated_day7_en_scores.csv"), show_col_types = FALSE)
chk <- calls %>% inner_join(old, by = "Core_Sample_Id") %>%
  mutate(same = (Call == "Predicted non-survivor") ==
                (Predicted_Group_medianA == "Predicted_NonSurvivor"))
stopifnot(nrow(chk) == 11, all(chk$same))
message("Split check: arm-centered median split == Predicted_Group_medianA for all 11 mice.")

## ---- the three contrasts ----
reb <- read_csv(file.path(REB, "results", "gsea_all_contrasts_long.csv"), show_col_types = FALSE)
pull_reb <- function(cmp) reb %>% filter(comparison == cmp) %>% dplyr::select(pathway, NES, padj)
pred <- read_csv(file.path(SRC, "results", "gsea_day7_untreated_pred_median_long.csv"),
                 show_col_types = FALSE) %>% dplyr::select(pathway, NES, padj)
ut   <- pull_reb("Day7_untreated")
tr   <- pull_reb("Day7_treated")

## Sign convention is shared across all three: NES > 0 = up in (predicted)
## non-survivors, matching the published Figure 4/6 panels.
CONTRASTS <- list(
  `EN-predicted split\n(untreated, n=11)`  = pred,
  `Actual survival\n(untreated, n=11)`     = ut,
  `Actual survival\n(treated, n=12)`       = tr)

catalog <- og_full_set %>% mutate(module = factor(module, levels = CAT_LEVELS))
long <- bind_rows(lapply(names(CONTRASTS), function(nm)
  CONTRASTS[[nm]] %>% inner_join(catalog, by = "pathway") %>% mutate(contrast = nm)))
stopifnot(all(table(long$contrast) == nrow(catalog)))   # no blank cells
write_csv(long, file.path(OUT, "tables", "S_gsea_catalog_NES.csv"))

## ---- alignment: which real axis does the predicted split look like? ----
cosim <- function(a, b) sum(a * b) / sqrt(sum(a^2) * sum(b^2))
sgn   <- function(a, b) mean(sign(a) == sign(b))

wide_all <- pred %>% rename(NES_pred = NES, padj_pred = padj) %>%
  inner_join(ut %>% rename(NES_ut = NES, padj_ut = padj), by = "pathway") %>%
  inner_join(tr %>% rename(NES_tr = NES, padj_tr = padj), by = "pathway")
wide_cat <- wide_all %>% filter(pathway %in% catalog$pathway)

align <- bind_rows(
  data.frame(scope = "31-pathway manuscript catalog", n = nrow(wide_cat),
             comparison = "EN-predicted split vs ACTUAL survival, treated",
             cosine = cosim(wide_cat$NES_pred, wide_cat$NES_tr),
             sign_concordance = sgn(wide_cat$NES_pred, wide_cat$NES_tr)),
  data.frame(scope = "31-pathway manuscript catalog", n = nrow(wide_cat),
             comparison = "EN-predicted split vs ACTUAL survival, untreated",
             cosine = cosim(wide_cat$NES_pred, wide_cat$NES_ut),
             sign_concordance = sgn(wide_cat$NES_pred, wide_cat$NES_ut)),
  data.frame(scope = "31-pathway manuscript catalog", n = nrow(wide_cat),
             comparison = "reference: treated vs untreated ACTUAL survival",
             cosine = cosim(wide_cat$NES_tr, wide_cat$NES_ut),
             sign_concordance = sgn(wide_cat$NES_tr, wide_cat$NES_ut)),
  data.frame(scope = "all 1,187 gene sets", n = nrow(wide_all),
             comparison = "EN-predicted split vs ACTUAL survival, treated",
             cosine = cosim(wide_all$NES_pred, wide_all$NES_tr),
             sign_concordance = sgn(wide_all$NES_pred, wide_all$NES_tr)),
  data.frame(scope = "all 1,187 gene sets", n = nrow(wide_all),
             comparison = "EN-predicted split vs ACTUAL survival, untreated",
             cosine = cosim(wide_all$NES_pred, wide_all$NES_ut),
             sign_concordance = sgn(wide_all$NES_pred, wide_all$NES_ut)),
  data.frame(scope = "all 1,187 gene sets", n = nrow(wide_all),
             comparison = "reference: treated vs untreated ACTUAL survival",
             cosine = cosim(wide_all$NES_tr, wide_all$NES_ut),
             sign_concordance = sgn(wide_all$NES_tr, wide_all$NES_ut)))
print(align, digits = 3)
write_csv(align, file.path(OUT, "tables", "S_gsea_axis_alignment.csv"))

## ---- variant sensitivity: does the conclusion depend on how the score was cut? ----
## The three alternative splits from 02_gsea_predicted_groups.R, scored the same way.
variants <- c(`Median split (5v6, n=11)`      = "gsea_day7_untreated_pred_median_long.csv",
              `Natural-gap split (9v2, n=11)` = "gsea_day7_untreated_pred_gap_long.csv",
              `Margin-excluded (4v5, n=9)`    = "gsea_day7_untreated_pred_margin_long.csv",
              `Continuous score (n=11)`       = "gsea_day7_untreated_pred_continuous_long.csv")
sens <- bind_rows(lapply(names(variants), function(v) {
  x <- read_csv(file.path(SRC, "results", variants[[v]]), show_col_types = FALSE) %>%
    dplyr::select(pathway, NES) %>%
    inner_join(tr %>% rename(NES_tr = NES), by = "pathway") %>%
    inner_join(ut %>% rename(NES_ut = NES), by = "pathway")
  xc <- x %>% filter(pathway %in% catalog$pathway)
  data.frame(variant = v,
             cat31_cos_vs_treated   = cosim(xc$NES, xc$NES_tr),
             cat31_cos_vs_untreated = cosim(xc$NES, xc$NES_ut),
             all_cos_vs_treated     = cosim(x$NES,  x$NES_tr),
             all_cos_vs_untreated   = cosim(x$NES,  x$NES_ut))
}))
print(sens, digits = 3)
write_csv(sens, file.path(OUT, "tables", "S_gsea_split_sensitivity.csv"))
message("Wrote GSEA tables.")
