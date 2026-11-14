## Panel helpers for the untreated-mouse supplementary figure (heatmaps, ROC,
## risk scores, permutation nulls).

suppressPackageStartupMessages({
  .libPaths(c("./NCT_paper_pipeline/renv_lib",
              .libPaths()))
  library(readr); library(dplyr); library(matrixStats); library(limma)
  library(glmnet); library(pROC); library(pheatmap); library(ggplot2)
  library(edgeR); library(grid); library(gridExtra); library(png)
})

set.seed(123)

PROJ        <- "."
PIPE        <- file.path(PROJ, "NCT_paper_pipeline")
INPUT_COUNTS<- file.path(PIPE, "01_inputs", "merged_counts.tsv")
INPUT_META  <- file.path(PIPE, "01_inputs", "merged_metadata_all.tsv")
D18_ROOT <- file.path(PIPE, "02_cache", "elastic_net_day7_day18_source_match",
                      "results_EN_day18_untreated_batch1_contam_filtered_15gene")
D7_ROOT  <- file.path(PIPE, "02_cache", "elastic_net_day7_day18_source_match",
                      "results_EN_day7_untreated_contam_filtered_15gene")
REPORT_ROOT <- file.path(PROJ, "Figure3_signature_analysis_report")
FIG_DIR     <- "./Figure3_signature_analysis_report/figures/suppfig_untreated"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

COLORS_BIN <- c(Survivor = "#1f77b4", NonSurvivor = "#d62728")
COLORS_TP  <- c("Day 18 (resection)" = "#117733", "Day 7 (pre-treatment)" = "#882255")

## ---- helpers (identical to 09/11/12/15) ----
strip_ver  <- function(x) sub("\\.\\d+$", "", x)
zscore_rows<- function(mat){ z <- t(scale(t(mat))); z[is.na(z)] <- 0; z }
map_sym <- function(ens) {
  ens <- strip_ver(ens)
  if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    odb <- getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
    tab <- AnnotationDbi::select(odb, keys = ens, keytype = "ENSEMBL", columns = "SYMBOL")
    out <- rep(NA_character_, length(ens))
    out[match(tab$ENSEMBL, ens)] <- tab$SYMBOL
    return(out)
  }
  rep(NA_character_, length(ens))
}

meta_all <- read_tsv(INPUT_META, show_col_types = FALSE)
cts_tbl  <- read_tsv(INPUT_COUNTS, show_col_types = FALSE)
gene_col <- names(cts_tbl)[1]
sample_cols <- setdiff(names(cts_tbl), gene_col)
ens <- strip_ver(as.character(cts_tbl[[gene_col]]))
mat_raw <- as.matrix(cts_tbl[, sample_cols]); storage.mode(mat_raw) <- "numeric"
rownames(mat_raw) <- ens
mat_raw <- rowsum(mat_raw, group = rownames(mat_raw), reorder = FALSE)
sym <- map_sym(rownames(mat_raw))
keep_sym <- !is.na(sym) & nzchar(sym) & !grepl("Rik$", sym) & !grepl("^(Gm|Ps)", sym)
mat_raw <- mat_raw[keep_sym, ]
mat_sym <- rowsum(mat_raw, group = sym[keep_sym], reorder = FALSE)

build_inputs <- function(timepoint, treatment_filter, do_batch_remove, batch_filter = NULL) {
  m <- meta_all %>%
    filter(Timepoint == timepoint,
           Survival_Status %in% c("survivor", "non-survivor")) %>%
    mutate(Batch = factor(as.character(Batch)),
           Treatment = tolower(Treatment),
           GroupBin = ifelse(tolower(Survival_Status) == "non-survivor",
                             "NonSurvivor", "Survivor"))
  if (!is.null(treatment_filter)) m <- m %>% filter(Treatment == treatment_filter)
  if (!is.null(batch_filter)) m <- m %>% filter(as.character(Batch) == batch_filter)
  common <- intersect(colnames(mat_sym), m$Core_Sample_Id)
  Xraw   <- mat_sym[, common, drop = FALSE]
  m      <- m %>% filter(Core_Sample_Id %in% common) %>%
    arrange(match(Core_Sample_Id, common)) %>% as.data.frame()
  d <- DGEList(counts = Xraw)
  d <- calcNormFactors(d, method = "TMM")
  cpm_mat <- cpm(d, log = FALSE)
  keep <- rowSums(cpm_mat >= 1) >= ceiling(0.25 * ncol(cpm_mat))
  d <- d[keep, , keep.lib.sizes = FALSE]
  mat <- cpm(d, log = TRUE, prior.count = 1)
  if (do_batch_remove && nlevels(m$Batch) > 1) mat <- removeBatchEffect(mat, batch = m$Batch)
  rv <- rowVars(mat); mat <- mat[rv > 0, ]
  list(mat_z = zscore_rows(mat), meta = m,
       y = as.integer(m$GroupBin == "NonSurvivor"))
}

## ---- per-contrast assembly ----
assemble <- function(sig_root, contrast_dir, timepoint, do_batch_remove, tp_label,
                     treatment_filter = "n", batch_filter = NULL) {
  inp  <- build_inputs(timepoint, treatment_filter, do_batch_remove, batch_filter)
  best <- read_csv(file.path(sig_root, contrast_dir, "selected_best_candidate.csv"),
                   show_col_types = FALSE)
  sig_tab <- read_csv(file.path(sig_root, contrast_dir, "selected_signature_genes.csv"),
                      show_col_types = FALSE) %>% arrange(desc(abs(coef)))
  sig_genes <- intersect(sig_tab$gene, rownames(inp$mat_z))
  scores <- read_csv(file.path(sig_root, contrast_dir, "selected_LOOCV_scores.csv"),
                     show_col_types = FALSE)
  scores <- scores[match(inp$meta$Core_Sample_Id, scores$Core_Sample_Id), ]
  roc_obj <- suppressWarnings(roc(inp$y, scores$LOOCV_score,
                                  levels = c(0,1), direction = "<", quiet = TRUE))
  ci <- suppressWarnings(ci.auc(roc_obj, method = "delong"))
  perm <- read_csv(file.path(sig_root, contrast_dir, "validation_ci_perm",
                             "permutation_null_aucs.csv"), show_col_types = FALSE)
  vsum_path_top    <- file.path(sig_root, "validation_ci_perm_summary.csv")
  vsum_path_nested <- file.path(sig_root, contrast_dir, "validation_ci_perm",
                                "validation_ci_perm_summary.csv")
  if (file.exists(vsum_path_top)) {
    vsum <- read_csv(vsum_path_top, show_col_types = FALSE) %>% filter(Contrast == contrast_dir)
  } else {
    vsum <- read_csv(vsum_path_nested, show_col_types = FALSE)
  }
  # Refit-on-all (apparent) risk scores for panel D: apply the signature
  # coefficients to the z-scored matrix and calibrate the intercept by ML.
  sg_coef <- sig_tab$coef[match(sig_genes, sig_tab$gene)]
  eta <- as.numeric(t(inp$mat_z[sig_genes, , drop = FALSE]) %*% sg_coef)
  b0  <- suppressWarnings(coef(glm(inp$y ~ 1 + offset(eta), family = binomial))[1])
  refit_scores <- plogis(b0 + eta)
  list(meta = inp$meta, mat_z = inp$mat_z, sig_genes = sig_genes,
       scores = scores$LOOCV_score, refit_scores = refit_scores,
       group = inp$meta$GroupBin,
       roc = roc_obj, auc = as.numeric(auc(roc_obj)), ci = ci,
       null_auc = perm$null_auc, perm_p = vsum$Perm_p[1], n = length(inp$y),
       K = best$Selected_K[1], alpha = best$Selected_alpha[1],
       tp_label = tp_label)
}

A <- assemble(D18_ROOT, "Day18_Untreated_SvNS", "Day18", FALSE, "Day 18 (resection)",
             batch_filter = "batch_1")
B <- assemble(D7_ROOT, "Day7_Untreated_SvNS", "Day7", FALSE, "Day 7 (pre-treatment)")

message(sprintf("Panel A (Day18, batch_1, contam-filtered): N=%d, K=%d, alpha=%.2f, genes=%d, AUC=%.3f, perm p=%.3f",
                A$n, A$K, A$alpha, length(A$sig_genes), A$auc, A$perm_p))
message(sprintf("Panel B (Day7, contam-filtered):           N=%d, K=%d, alpha=%.2f, genes=%d, AUC=%.3f, perm p=%.3f",
                B$n, B$K, B$alpha, length(B$sig_genes), B$auc, B$perm_p))


## ==============================================================
## panel assembly - no PCA, Day 18 risk scores only
## ==============================================================
theme_fig <- theme_classic(base_size = 11) +
  theme(plot.margin = margin(4,6,4,6), legend.key.size = unit(0.8,"lines"))

render_heatmap_clean <- function(d, file) {
  M <- d$mat_z[d$sig_genes, , drop=FALSE]
  ann <- data.frame(Group = factor(d$group, levels=c("Survivor","NonSurvivor"))); rownames(ann) <- colnames(M)
  png(file, width=1350, height=1250, res=200)
  pheatmap(M, scale="none", show_colnames=FALSE, show_rownames=TRUE, fontsize_row=9, fontsize=9,
           clustering_method="complete", main="", annotation_col=ann,
           annotation_colors=list(Group=setNames(c(COLORS_BIN["Survivor"],COLORS_BIN["NonSurvivor"]),
                                                  c("Survivor","NonSurvivor"))),
           annotation_names_col=FALSE)
  dev.off()
}
heatmap_grob <- function(file) rasterGrob(readPNG(file), width=unit(1,"npc"), height=unit(1,"npc"))

roc_plot_clean <- function() {
  rd <- function(d, lab){ rc <- as.data.frame(coords(d$roc,"all",ret=c("specificity","sensitivity"),transpose=FALSE))
    out <- data.frame(FPR=1-rc$specificity, TPR=rc$sensitivity, Sig=lab); out[order(out$FPR,out$TPR),] }
  df <- rbind(rd(A,"Day 18 (resection)"), rd(B,"Day 7 (pre-treatment)"))
  cols <- c("Day 18 (resection)"="#117733","Day 7 (pre-treatment)"="#882255")
  ggplot(df, aes(FPR, TPR, color=Sig)) +
    geom_abline(linetype=2, color="grey60") + geom_path(linewidth=1.1) +
    scale_color_manual(values=cols, guide="none") + coord_equal() +
    annotate("text", x=0.97, y=0.16, hjust=1, label=sprintf("D18 AUC= %.2f", A$auc),
             color=cols[1], size=3.6) +
    annotate("text", x=0.97, y=0.06, hjust=1, label=sprintf("D7  AUC= %.2f", B$auc),
             color=cols[2], size=3.6) +
    labs(x="False positive rate", y="True positive rate") + theme_fig
}

## Day 18 only, with the Wilcoxon p printed
scores_plot_day18 <- function() {
  sc <- data.frame(Score=A$refit_scores, Group=factor(A$group, levels=c("Survivor","NonSurvivor")))
  pv <- suppressWarnings(wilcox.test(Score ~ Group, data=sc)$p.value)
  ggplot(sc, aes(Group, Score, fill=Group)) +
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.6) +
    geom_jitter(aes(color=Group), width=0.12, size=2.1, alpha=0.9) +
    annotate("text", x=1.5, y=1.06,
             label=if (pv < 0.001) "Wilcoxon p < 0.001" else sprintf("Wilcoxon p = %.3f", pv),
             size=3.4) +
    scale_fill_manual(values=COLORS_BIN[c("Survivor","NonSurvivor")]) +
    scale_color_manual(values=COLORS_BIN[c("Survivor","NonSurvivor")]) +
    coord_cartesian(ylim=c(0,1.12)) +
    labs(x=NULL, y="EN score (P[non-survival])") +
    theme_fig + theme(legend.position="none")
}

null_plot_clean <- function(d, tag, lab) {
  ggplot(data.frame(x=d$null_auc), aes(x)) +
    geom_histogram(bins=25, fill="grey78", color="white") +
    geom_vline(xintercept=d$auc, color=COLORS_BIN["NonSurvivor"], linewidth=1.1) +
    annotate("text", x=d$auc, y=Inf, vjust=1.6, hjust=1.05,
             label=sprintf("AUC %.2f
p= %.3f", d$auc, d$perm_p), size=3, color=COLORS_BIN["NonSurvivor"]) +
    coord_cartesian(xlim=c(0,1)) +
    labs(tag=tag, x="Null AUC (permuted labels)", y="Count", subtitle=lab) +
    theme_fig + theme(plot.tag=element_text(size=12), plot.tag.position=c(0.01,0.99),
                      plot.subtitle=element_text(size=9, hjust=0.5))
}

hdr <- function(letter, title) gTree(children = gList(
  textGrob(letter, x=unit(0.5,"lines"), y=0.5, just=c("left","center"), gp=gpar(fontsize=16, fontface="bold")),
  textGrob(title,  x=0.5, y=0.5, just="center", gp=gpar(fontsize=12, fontface="bold"))))
HH <- unit(1.5, "lines")

graphics.off()
tmpf <- function(n) file.path(FIG_DIR, paste0("_hm_", n, ".png"))
invisible(render_heatmap_clean(A, tmpf("A")))
invisible(render_heatmap_clean(B, tmpf("B")))

panelA <- arrangeGrob(hdr("A", "Resection day (untreated, batch 1)"),
                      heatmap_grob(tmpf("A")), heights = unit.c(HH, unit(1,"null")))
panelB <- arrangeGrob(hdr("B", "Pre-treatment (untreated)"),
                      heatmap_grob(tmpf("B")), heights = unit.c(HH, unit(1,"null")))
panelC <- arrangeGrob(hdr("C", "Cross-validated ROC"), roc_plot_clean(),
                      heights=unit.c(HH, unit(1,"null")))
panelD <- arrangeGrob(hdr("D", "Risk scores, resection day"), scores_plot_day18(),
                      heights=unit.c(HH, unit(1,"null")))
panelE <- arrangeGrob(hdr("E", "Permutation validation"),
                      arrangeGrob(null_plot_clean(A,"i.)","Day 18 (resection), batch 1"),
                                  null_plot_clean(B,"ii.)","Day 7 (pre-treatment)"), ncol=2),
                      heights = unit.c(HH, unit(1,"null")))
rowAB <- arrangeGrob(panelA, panelB, ncol=2)
rowCD <- arrangeGrob(panelC, panelD, ncol=2)

fig <- arrangeGrob(rowAB, rowCD, panelE, ncol=1, heights=c(1.20,1.05,0.90))
png(file.path(FIG_DIR, "SuppFig_S1_noPCA_Day18scores.png"), width=2400, height=2500, res=215)
grid.draw(fig); dev.off()
file.remove(list.files(FIG_DIR, pattern="^_hm_.*png$", full.names=TRUE))

cat("
Wrote:", file.path(FIG_DIR, "SuppFig_S1_noPCA_Day18scores.png"), "
")
cat(sprintf("Day 18: n=%d, %d genes, AUC %.3f (%.3f-%.3f), perm p=%.3f
",
            A$n, length(A$sig_genes), A$auc, A$ci[1], A$ci[3], A$perm_p))
cat(sprintf("Day  7: n=%d, %d genes, AUC %.3f (%.3f-%.3f), perm p=%.3f
",
            B$n, length(B$sig_genes), B$auc, B$ci[1], B$ci[3], B$perm_p))
