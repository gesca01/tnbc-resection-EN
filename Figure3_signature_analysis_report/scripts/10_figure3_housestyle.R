## Figure 3 panels (treated mice): signature heatmaps, cross-validated ROC,
## risk scores and permutation nulls, from the cached signature run.

suppressPackageStartupMessages({
  .libPaths(c("./NCT_paper_pipeline/renv_lib",
              .libPaths()))
  library(readr); library(dplyr); library(matrixStats); library(limma)
  library(glmnet); library(pROC); library(pheatmap); library(ggplot2)
  library(edgeR); library(grid); library(gridExtra); library(png)
})
set.seed(123)

PROJ <- "."; PIPE <- file.path(PROJ, "NCT_paper_pipeline")
INPUT_COUNTS <- file.path(PIPE, "01_inputs", "merged_counts.tsv")
INPUT_META   <- file.path(PIPE, "01_inputs", "merged_metadata_all.tsv")
SIG_ROOT <- file.path(PIPE, "02_cache", "elastic_net_day7_day18_source_match",
                      "results_EN_signatures")
FIG_DIR <- file.path(PROJ, "Figure3_signature_analysis_report", "figures", "manuscript_treated")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
COLORS_BIN <- c(Survivor = "#1f77b4", NonSurvivor = "#d62728")

## ---- preprocessing (same as run_day7_day18_signatures.R) ----
strip_ver   <- function(x) sub("\\.\\d+$", "", x)
zscore_rows <- function(m){ z <- t(scale(t(m))); z[is.na(z)] <- 0; z }
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

build_inputs <- function(timepoint, do_batch_remove) {
  m <- meta_all %>% filter(Timepoint == timepoint, Survival_Status %in% c("survivor","non-survivor")) %>%
    mutate(Batch=factor(as.character(Batch)), Treatment=tolower(Treatment),
           GroupBin=ifelse(tolower(Survival_Status)=="non-survivor","NonSurvivor","Survivor")) %>%
    filter(Treatment == "y")
  common <- intersect(colnames(mat_sym), m$Core_Sample_Id); Xraw <- mat_sym[, common, drop=FALSE]
  m <- m %>% filter(Core_Sample_Id %in% common) %>% arrange(match(Core_Sample_Id, common)) %>% as.data.frame()
  d <- DGEList(Xraw); d <- calcNormFactors(d, "TMM"); cpm_mat <- cpm(d, log=FALSE)
  keep <- rowSums(cpm_mat >= 1) >= ceiling(0.25*ncol(cpm_mat)); d <- d[keep,,keep.lib.sizes=FALSE]
  mat <- cpm(d, log=TRUE, prior.count=1)
  if (do_batch_remove && nlevels(m$Batch) > 1) mat <- removeBatchEffect(mat, batch=m$Batch)
  mat <- mat[rowVars(mat) > 0, ]; list(mat_z=zscore_rows(mat), meta=m, y=as.integer(m$GroupBin=="NonSurvivor"))
}
assemble <- function(contrast_dir, timepoint, do_batch_remove) {
  inp <- build_inputs(timepoint, do_batch_remove)
  sig_tab <- read_csv(file.path(SIG_ROOT, contrast_dir, "selected_signature_genes.csv"),
                  show_col_types=FALSE) %>% arrange(desc(abs(coef)))
  g <- intersect(sig_tab$gene, rownames(inp$mat_z)); b <- sig_tab$coef[match(g, sig_tab$gene)]
  sc <- read_csv(file.path(SIG_ROOT, contrast_dir, "selected_LOOCV_scores.csv"), show_col_types=FALSE)
  sc <- sc[match(inp$meta$Core_Sample_Id, sc$Core_Sample_Id), ]
  roc_obj <- suppressWarnings(roc(inp$y, sc$LOOCV_score, levels=c(0,1), direction="<", quiet=TRUE))
  ci <- suppressWarnings(ci.auc(roc_obj, method="delong"))
  perm <- read_csv(file.path(SIG_ROOT, contrast_dir, "validation_ci_perm","permutation_null_aucs.csv"),
                   show_col_types=FALSE)
  vsum <- read_csv(file.path(SIG_ROOT,"validation_ci_perm_summary.csv"), show_col_types=FALSE) %>%
    filter(Contrast==contrast_dir)
  eta <- as.numeric(t(inp$mat_z[g,,drop=FALSE]) %*% b)
  b0  <- suppressWarnings(coef(glm(inp$y ~ 1 + offset(eta), family=binomial))[1])
  list(meta=inp$meta, mat_z=inp$mat_z, sig=g, group=inp$meta$GroupBin,
       roc=roc_obj, auc=as.numeric(auc(roc_obj)), ci=ci,
       refit=plogis(b0+eta), null=perm$null_auc, perm_p=vsum$Perm_p[1], n=length(inp$y))
}
A <- assemble("Day18_Treated_SvNS","Day18",TRUE)   # 3A resection
B <- assemble("Day7_Treated_SvNS","Day7",FALSE)    # 3B pre-treatment

## ---- clean panel builders (no titles) ----
theme_fig <- theme_classic(base_size = 11) +
  theme(plot.margin = margin(4,6,4,6), legend.key.size = unit(0.8,"lines"))

# Render a clean heatmap PNG to disk (done up front, before any ggplot grob
# conversion, so the pheatmap device write is not disturbed).
render_heatmap <- function(d, file) {
  M <- d$mat_z[d$sig, , drop=FALSE]
  ann <- data.frame(Group = factor(d$group, levels=c("Survivor","NonSurvivor"))); rownames(ann) <- colnames(M)
  png(file, width=1450, height=1250, res=200)
  pheatmap(M, scale="none", show_colnames=FALSE, show_rownames=TRUE, fontsize_row=8, fontsize=9,
           clustering_method="complete", main="", annotation_col=ann,
           annotation_colors=list(Group=setNames(c(COLORS_BIN["Survivor"],COLORS_BIN["NonSurvivor"]),
                                                  c("Survivor","NonSurvivor"))),
           annotation_names_col=FALSE)
  dev.off()
}
heatmap_grob <- function(file, tag) gTree(children = gList(
  rasterGrob(readPNG(file), width=unit(1,"npc"), height=unit(1,"npc")),
  textGrob(tag, x=unit(0.4,"lines"), y=unit(1,"npc")-unit(0.4,"lines"),
           just=c("left","top"), gp=gpar(fontsize=12))))
pca_plot <- function(d, tag) {
  pcs <- prcomp(t(d$mat_z[d$sig,,drop=FALSE]), scale.=FALSE)
  v <- 100*(pcs$sdev^2/sum(pcs$sdev^2))
  df <- data.frame(PC1=pcs$x[,1], PC2=pcs$x[,2], Group=factor(d$group, levels=c("Survivor","NonSurvivor")))
  ggplot(df, aes(PC1, PC2, color=Group)) +
    geom_point(size=2.6, alpha=0.9) +
    stat_ellipse(aes(group=Group), level=0.68, linewidth=0.6, alpha=0.6) +
    scale_color_manual(values=COLORS_BIN[c("Survivor","NonSurvivor")]) +
    labs(tag=tag, x=sprintf("PC1 (%.0f%%)", v[1]), y=sprintf("PC2 (%.0f%%)", v[2])) +
    theme_fig + theme(legend.position="top", legend.title=element_blank(),
                      plot.tag=element_text(size=12), plot.tag.position=c(0.01,0.99))
}
roc_plot <- function() {
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
scores_plot <- function() {
  sc <- rbind(data.frame(TP="Day 18 (resection)", Score=A$refit, Group=A$group),
              data.frame(TP="Day 7 (pre-treatment)", Score=B$refit, Group=B$group))
  sc$TP <- factor(sc$TP, levels=c("Day 7 (pre-treatment)","Day 18 (resection)"))
  sc$Group <- factor(sc$Group, levels=c("Survivor","NonSurvivor"))
  ggplot(sc, aes(Group, Score, fill=Group)) +
    geom_boxplot(width=0.55, outlier.shape=NA, alpha=0.6) +
    geom_jitter(aes(color=Group), width=0.12, size=1.8, alpha=0.9) +
    facet_wrap(~TP) +
    scale_fill_manual(values=COLORS_BIN[c("Survivor","NonSurvivor")]) +
    scale_color_manual(values=COLORS_BIN[c("Survivor","NonSurvivor")]) +
    labs(x=NULL, y="EN score (P[non-survival])") +
    theme_fig + theme(legend.position="none")
}
null_plot <- function(d, tag, lab) {
  ggplot(data.frame(x=d$null), aes(x)) +
    geom_histogram(bins=25, fill="grey78", color="white") +
    geom_vline(xintercept=d$auc, color=COLORS_BIN["NonSurvivor"], linewidth=1.1) +
    annotate("text", x=d$auc, y=Inf, vjust=1.6, hjust=1.05,
             label=sprintf("AUC %.2f\np= %.3f", d$auc, d$perm_p), size=3, color=COLORS_BIN["NonSurvivor"]) +
    coord_cartesian(xlim=c(0,1)) +
    labs(tag=tag, x="Null AUC (permuted labels)", y="Count", subtitle=lab) +
    theme_fig + theme(plot.tag=element_text(size=12), plot.tag.position=c(0.01,0.99),
                      plot.subtitle=element_text(size=9, hjust=0.5))
}

## ---- house-style header (bold letter at left, title centered) ----
hdr <- function(letter, title) gTree(children = gList(
  textGrob(letter, x=unit(0.5,"lines"), y=0.5, just=c("left","center"), gp=gpar(fontsize=16, fontface="bold")),
  textGrob(title,  x=0.5, y=0.5, just="center", gp=gpar(fontsize=12, fontface="bold"))))
HH <- unit(1.5, "lines")   # header height

tmpf <- function(n) file.path(FIG_DIR, paste0("_hm_", n, ".png"))
render_heatmap(A, tmpf("A")); render_heatmap(B, tmpf("B"))   # render up front
panelA <- arrangeGrob(hdr("A", "Resection-day signature (treated)"),
                      arrangeGrob(heatmap_grob(tmpf("A"),"i.)"), pca_plot(A,"ii.)"), ncol=2),
                      heights = unit.c(HH, unit(1,"null")))
panelB <- arrangeGrob(hdr("B", "Pre-treatment signature (treated)"),
                      arrangeGrob(heatmap_grob(tmpf("B"),"i.)"), pca_plot(B,"ii.)"), ncol=2),
                      heights = unit.c(HH, unit(1,"null")))
panelC <- arrangeGrob(hdr("C", "Cross-validated ROC"), roc_plot(), heights=unit.c(HH, unit(1,"null")))
panelD <- arrangeGrob(hdr("D", "Risk scores by group"), scores_plot(), heights=unit.c(HH, unit(1,"null")))
panelE <- arrangeGrob(hdr("E", "Permutation validation"),
                      arrangeGrob(null_plot(A,"i.)","Day 18 (resection)"),
                                  null_plot(B,"ii.)","Day 7 (pre-treatment)"), ncol=2),
                      heights = unit.c(HH, unit(1,"null")))
rowCD <- arrangeGrob(panelC, panelD, ncol=2)

fig <- arrangeGrob(panelA, panelB, rowCD, panelE, ncol=1, heights=c(1.18,1.18,1.05,0.95))
png(file.path(FIG_DIR, "Figure3_treated_housestyle.png"), width=2500, height=4150, res=215)
grid.draw(fig); dev.off()
file.remove(list.files(FIG_DIR, pattern="^_hm_.*png$", full.names=TRUE))
message("Done: ", file.path(FIG_DIR, "Figure3_treated_housestyle.png"))
