## GSEA on untreated day-7 mice grouped by the outcome the day-7 treated signature
## predicts for them (from 01), instead of their observed outcome.
## Same pipeline as run_gsea_all_contrasts.R (DESeq2 Wald stat -> fgsea).
## Contrasts: predicted groups split at the median score (5 vs 6), split at the
## largest gap in scores (9 vs 2), and observed outcome as a check.
## NES > 0 = higher in (predicted) non-survivors.

options(stringsAsFactors = FALSE)
set.seed(123)

PROJ <- "."
PIPELINE_ROOT <- file.path(PROJ, "NCT_paper_pipeline")
OUT <- file.path(PROJ, "Day7_Treated_Signature_Applied_to_Untreated")
RES <- file.path(OUT, "results")
.libPaths(file.path(PIPELINE_ROOT, "renv_lib"))

COUNTS_TSV   <- file.path(PIPELINE_ROOT, "01_inputs", "merged_counts.tsv")
METADATA_TSV <- file.path(PIPELINE_ROOT, "01_inputs", "merged_metadata_all.tsv")
SCORES_CSV   <- file.path(RES, "untreated_day7_en_scores.csv")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr); library(tibble)
  library(DESeq2); library(fgsea); library(msigdbr); library(AnnotationDbi); library(org.Mm.eg.db)
})
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

## ---- helpers (same as run_gsea_all_contrasts.R) ----
strip_version <- function(x) sub("\\.\\d+$", "", x)
mk_rank <- function(res) { r <- res$stat; names(r) <- rownames(res); sort(r[is.finite(r)], decreasing = TRUE) }
run_fgsea <- function(ranks, pathways) {
  res <- fgsea(pathways, ranks, minSize = 15, maxSize = 500)
  as_tibble(res) %>% mutate(leadingEdgePct = round(100*lengths(leadingEdge)/size,1),
                            leadingEdge = vapply(leadingEdge, function(x) paste(unlist(x), collapse=","), character(1))) %>%
    arrange(padj, desc(abs(NES))) %>% select(pathway, size, ES, NES, pval, padj, leadingEdgePct, leadingEdge)
}
get_coef_name <- function(dds, pattern) { nm <- resultsNames(dds); hit <- nm[grepl(pattern, nm)]
  if (length(hit) != 1) stop("coef pattern ", pattern, " -> ", paste(nm, collapse=" | ")); hit }

## ---- load + recode exactly as the original does ----
message("Loading data ...")
counts_raw <- read_tsv(COUNTS_TSV, col_types = cols()); meta_raw <- read_tsv(METADATA_TSV, col_types = cols())
gene_col <- names(counts_raw)[1]; counts_raw <- counts_raw %>% rename(gene_id = all_of(gene_col))
meta <- meta_raw %>% mutate(across(c(Timepoint,Treatment,Survival_Status,Batch), as.character)) %>%
  filter(Survival_Status %in% c("survivor","non-survivor")) %>%
  mutate(Treatment = case_when(Treatment=="y"~"treated", Treatment=="n"~"untreated", TRUE~NA_character_),
         Survival_Status = case_when(Survival_Status=="non-survivor"~"non_survivor", Survival_Status=="survivor"~"survivor", TRUE~NA_character_),
         Batch = str_replace_all(Batch, "[^A-Za-z0-9_]", "_")) %>%
  mutate(Treatment=factor(Treatment,levels=c("untreated","treated")),
         Survival_Status=factor(Survival_Status,levels=c("survivor","non_survivor")), Batch=factor(Batch))
sample_cols <- intersect(names(counts_raw), meta$Core_Sample_Id)
meta <- meta %>% filter(Core_Sample_Id %in% sample_cols) %>%
  mutate(Core_Sample_Id=factor(Core_Sample_Id,levels=sample_cols)) %>% arrange(Core_Sample_Id)
meta$Core_Sample_Id <- as.character(meta$Core_Sample_Id)
count_df <- counts_raw %>% select(gene_id, all_of(meta$Core_Sample_Id))

message("Mapping genes ...")
ens0 <- strip_version(count_df$gene_id)
sym <- mapIds(org.Mm.eg.db, keys=ens0, keytype="ENSEMBL", column="SYMBOL", multiVals="first")
gk <- ifelse(!is.na(sym) & sym!="", as.character(sym), ens0)
count_df2 <- count_df %>% mutate(gene_key=gk) %>% select(gene_key, all_of(meta$Core_Sample_Id)) %>%
  group_by(gene_key) %>% summarise(across(everything(), ~sum(as.numeric(.), na.rm=TRUE)), .groups="drop")
cts <- as.matrix(count_df2[,-1]); rownames(cts) <- count_df2$gene_key
cts_f <- cts[rowSums(cts >= 10) >= 3, , drop = FALSE]   # global filter, all samples
message("Genes kept: ", nrow(cts_f))

message("Building gene sets (H + Reactome + KEGG legacy) ...")
gs <- function(cat, sub=NULL) { a <- list(species="Mus musculus", category=cat); if(!is.null(sub)) a$subcategory <- sub
  do.call(msigdbr, a) %>% select(gs_name, gene_symbol) }
msig <- bind_rows(
  gs("H") %>% mutate(collection="HALLMARK"),
  gs("C2","CP:REACTOME") %>% mutate(collection="REACTOME"),
  gs("C2","CP:KEGG") %>% mutate(collection="KEGG"))
set_collection <- msig %>% distinct(gs_name, collection)
GDF <- msig %>% group_by(gs_name) %>% summarise(genes=list(unique(gene_symbol)), .groups="drop") %>%
  mutate(size=lengths(genes)) %>% filter(size>=15, size<=500)
GENESETS <- setNames(GDF$genes, GDF$gs_name)
message("Pathways: ", length(GENESETS))

## ---- untreated Day 7 subset + the EN-predicted grouping ----
m <- meta %>% filter(Timepoint == "Day7", Treatment == "untreated")
scores <- read_csv(SCORES_CSV, show_col_types = FALSE)
stopifnot(setequal(m$Core_Sample_Id, scores$Core_Sample_Id))
i <- match(m$Core_Sample_Id, scores$Core_Sample_Id)
lv <- c("Predicted_Survivor", "Predicted_NonSurvivor")
m$Pred_Median <- factor(scores$Predicted_Group_medianA[i], levels = lv)  # 5 vs 6, n=11
m$Pred_Gap    <- factor(scores$Predicted_Group_p50A[i],    levels = lv)  # 9 vs 2, n=11
m$Pred_Margin <- factor(scores$Predicted_Group_margin[i],  levels = lv)  # 4 vs 5, n=9 (NA = dropped)
m$eta_z       <- scores$eta_A_z[i]                                       # continuous, no split
m <- m %>% mutate(Batch = droplevels(Batch), Survival_Status = droplevels(Survival_Status))
cnt <- cts_f[, m$Core_Sample_Id, drop = FALSE]
message("Untreated Day7: n=", nrow(m),
        "  median ", sum(m$Pred_Median==lv[1]), "/", sum(m$Pred_Median==lv[2]),
        "  gap ", sum(m$Pred_Gap==lv[1]), "/", sum(m$Pred_Gap==lv[2]),
        "  margin ", sum(m$Pred_Margin==lv[1], na.rm=TRUE), "/", sum(m$Pred_Margin==lv[2], na.rm=TRUE),
        " (", sum(is.na(m$Pred_Margin)), " dropped)",
        "  actual ", sum(m$Survival_Status=="survivor"), "/", sum(m$Survival_Status=="non_survivor"),
        "  batch levels=", nlevels(m$Batch))

## Batch is constant here (all Day7 is batch_2), so it is correctly absent from
## the design - same reduction the original's build_design() performs.
stopifnot(nlevels(m$Batch) == 1)

## `continuous = TRUE` puts the score in as a numeric covariate - DESeq2 then
## reports the Wald statistic per 1-SD change in score, with no threshold
## anywhere and every mouse contributing. Samples with NA in the grouping
## variable (margin exclusion) are dropped before fitting.
run_contrast <- function(var, tag, desc, continuous = FALSE) {
  keep <- !is.na(m[[var]])
  md <- droplevels(m[keep, , drop = FALSE]); cd <- cnt[, keep, drop = FALSE]
  design <- as.formula(paste("~", var))
  dds <- DESeq(DESeqDataSetFromMatrix(round(cd), md, design))
  message("  resultsNames: ", paste(resultsNames(dds), collapse = " | "))
  res <- results(dds, name = get_coef_name(dds, if (continuous) paste0("^", var, "$")
                                                else paste0("^", var, "_")))
  g <- run_fgsea(mk_rank(res), GENESETS) %>%
    left_join(set_collection, by = c("pathway" = "gs_name")) %>%
    mutate(comparison = tag)
  message("Done ", tag, "  n=", nrow(md), "  sig05=", sum(g$padj < 0.05, na.rm = TRUE))
  list(long = g,
       manifest = tibble(tag = tag, n = nrow(md), design = deparse(design),
                         n_group1 = if (continuous) NA_integer_ else sum(md[[var]] == levels(md[[var]])[1]),
                         n_group2 = if (continuous) NA_integer_ else sum(md[[var]] == levels(md[[var]])[2]),
                         n_sig05 = sum(g$padj < 0.05, na.rm = TRUE), desc = desc))
}

set.seed(123)
p_med  <- run_contrast("Pred_Median", "Day7_untreated_pred_median",
                       "EN-predicted groups, median/rank split (5 vs 6, n=11)")
p_gap  <- run_contrast("Pred_Gap", "Day7_untreated_pred_gap",
                       "EN-predicted groups, natural-gap / p>0.5 split (9 vs 2, n=11)")
p_marg <- run_contrast("Pred_Margin", "Day7_untreated_pred_margin",
                       "EN-predicted groups, median split with the 2 borderline mice dropped (4 vs 5, n=9)")
p_cont <- run_contrast("eta_z", "Day7_untreated_pred_continuous",
                       "EN score as a continuous covariate -- no split, all 11 mice",
                       continuous = TRUE)
act    <- run_contrast("Survival_Status", "Day7_untreated_actual",
                       "Observed survival outcome (verification vs gsea_rebuild)")

write_csv(p_med$long,  file.path(RES, "gsea_day7_untreated_pred_median_long.csv"))
write_csv(p_gap$long,  file.path(RES, "gsea_day7_untreated_pred_gap_long.csv"))
write_csv(p_marg$long, file.path(RES, "gsea_day7_untreated_pred_margin_long.csv"))
write_csv(p_cont$long, file.path(RES, "gsea_day7_untreated_pred_continuous_long.csv"))
write_csv(act$long,    file.path(RES, "gsea_day7_untreated_actual_verification_long.csv"))
write_csv(bind_rows(p_med$manifest, p_gap$manifest, p_marg$manifest, p_cont$manifest, act$manifest),
          file.path(RES, "gsea_contrast_manifest.csv"))

## How much do the candidate groupings differ, and how does each relate to the
## observed-outcome contrast? The last column is the one that matters: how much
## real survival biology each way of using the score recovers.
nes_cor <- function(a, b) {
  j <- inner_join(select(a, pathway, x = NES), select(b, pathway, y = NES), by = "pathway")
  cor(j$x, j$y)
}
variants <- list(median = p_med, gap = p_gap, margin = p_marg, continuous = p_cont)
agree <- tibble(
  variant = names(variants),
  n = vapply(variants, function(v) v$manifest$n, numeric(1)),
  n_sig05 = vapply(variants, function(v) v$manifest$n_sig05, numeric(1)),
  NES_cor_vs_actual = vapply(variants, function(v) round(nes_cor(v$long, act$long), 4), numeric(1)),
  NES_cor_vs_median = vapply(variants, function(v) round(nes_cor(v$long, p_med$long), 4), numeric(1)))
write_csv(agree, file.path(RES, "split_agreement_NES_correlations.csv"))
message("\n=== How each grouping compares ==="); print(as.data.frame(agree), row.names = FALSE)

## ---- verification: does the ported engine reproduce the existing result? ----
ref <- read_csv(file.path(PROJ, "Manuscript_Revision_2026-07", "gsea_rebuild", "results",
                          "gsea_all_contrasts_long.csv"), show_col_types = FALSE) %>%
  filter(comparison == "Day7_untreated") %>% select(pathway, NES_ref = NES, padj_ref = padj)
chk <- act$long %>% select(pathway, NES, padj) %>% inner_join(ref, by = "pathway")
message(sprintf("\nVERIFICATION vs gsea_rebuild Day7_untreated: %d shared pathways, NES r = %.5f, max |dNES| = %.4f",
                nrow(chk), cor(chk$NES, chk$NES_ref), max(abs(chk$NES - chk$NES_ref))))

message("\n=== MANIFEST ===")
print(as.data.frame(bind_rows(p_med$manifest, p_gap$manifest, p_marg$manifest,
                              p_cont$manifest, act$manifest))[
  , c("tag","n","design","n_group1","n_group2","n_sig05")], row.names = FALSE)
message("Written to ", RES)
