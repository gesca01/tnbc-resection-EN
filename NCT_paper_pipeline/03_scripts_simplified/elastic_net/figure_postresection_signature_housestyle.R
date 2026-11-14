## Figure 5 panel helpers: post-resection signature heatmap, PCA, ROC and
## EN-score trajectory, from the cache written by run_postresection_signature_v2.R.

suppressPackageStartupMessages({
  .libPaths(c("./NCT_paper_pipeline/renv_lib", .libPaths()))
  library(readr); library(dplyr); library(ggplot2)
  library(pheatmap); library(grid); library(gridExtra); library(png); library(pROC)
})
set.seed(123)

PROJ <- "."; PIPE <- file.path(PROJ, "NCT_paper_pipeline")
SUF     <- Sys.getenv("OUT_SUFFIX", "")
CACHE   <- file.path(PIPE, "02_cache", paste0("bmes_postresection_v2_housestyle", SUF))
FIG_DIR <- file.path(PIPE, "04_outputs", "quick_figures", paste0("figure5_housestyle", SUF))
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
COLORS_BIN <- c(Survivor = "#1f77b4", NonSurvivor = "#d62728")

## ---- load cache ----
inp     <- readRDS(file.path(CACHE, "figure_inputs.rds"))
summ    <- read_csv(file.path(CACHE, "summary.csv"), show_col_types = FALSE)
loocv   <- read_csv(file.path(CACHE, "selected_LOOCV_scores.csv"), show_col_types = FALSE)
traj    <- read_csv(file.path(CACHE, "trajectory_scores.csv"), show_col_types = FALSE)
trajp   <- read_csv(file.path(CACHE, "trajectory_pvals.csv"), show_col_types = FALSE)
nulls   <- read_csv(file.path(CACHE, "permutation_null_aucs.csv"), show_col_types = FALSE)

auc_sel <- summ$AUC_LOOCV[1]; perm_p <- summ$Perm_p[1]; n_genes <- summ$Signature_genes[1]
mat_z_sig <- inp$mat_z_sig
group <- factor(inp$meta_m$GroupBin, levels = c("Survivor", "NonSurvivor"))
phase <- factor(as.character(inp$meta_m$Phase), levels = c("Week1", "Mid", "Late"))
treat <- factor(ifelse(inp$meta_m$Treatment == "y", "treated", "untreated"), levels = c("treated", "untreated"))

theme_fig <- theme_classic(base_size = 11) +
  theme(plot.margin = margin(4, 6, 4, 6), legend.key.size = unit(0.8, "lines"))

## ---- heatmap (rendered to disk first) ----
render_heatmap <- function(file) {
  M <- mat_z_sig
  ann <- data.frame(Group = group, Phase = phase); rownames(ann) <- colnames(M)
  ann_cols <- list(Group = setNames(c(COLORS_BIN["Survivor"], COLORS_BIN["NonSurvivor"]), c("Survivor", "NonSurvivor")),
                   Phase = setNames(c("#66c2a5", "#fc8d62", "#8da0cb"), c("Week1", "Mid", "Late")))
  png(file, width = 1500, height = 1300, res = 200)
  pheatmap(M, scale = "none", show_colnames = FALSE, show_rownames = TRUE, fontsize_row = 7.5, fontsize = 9,
           clustering_method = "complete", main = "", annotation_col = ann,
           annotation_colors = ann_cols, annotation_names_col = FALSE)
  dev.off()
}
heatmap_grob <- function(file, tag) gTree(children = gList(
  rasterGrob(readPNG(file), width = unit(1, "npc"), height = unit(1, "npc")),
  textGrob(tag, x = unit(0.4, "lines"), y = unit(1, "npc") - unit(0.4, "lines"),
           just = c("left", "top"), gp = gpar(fontsize = 12))))

pca_plot <- function(tag) {
  pcs <- prcomp(t(mat_z_sig), scale. = FALSE); v <- 100 * (pcs$sdev^2 / sum(pcs$sdev^2))
  df <- data.frame(PC1 = pcs$x[, 1], PC2 = pcs$x[, 2], Group = group, Treatment = treat)
  ggplot(df, aes(PC1, PC2, color = Group, shape = Treatment)) +
    geom_point(size = 2.6, alpha = 0.9) +
    stat_ellipse(aes(group = Group), level = 0.68, linewidth = 0.6, alpha = 0.6) +
    scale_color_manual(values = COLORS_BIN[c("Survivor", "NonSurvivor")]) +
    labs(tag = tag, x = sprintf("PC1 (%.0f%%)", v[1]), y = sprintf("PC2 (%.0f%%)", v[2])) +
    theme_fig + theme(legend.position = "top", legend.title = element_blank(), legend.box = "horizontal",
                      plot.tag = element_text(size = 12), plot.tag.position = c(0.01, 0.99))
}
roc_plot <- function() {
  y <- as.integer(loocv$Group == "NonSurvivor")
  rocobj <- suppressWarnings(roc(y, loocv$LOOCV_score, levels = c(0, 1), direction = "<", quiet = TRUE))
  rc <- as.data.frame(coords(rocobj, "all", ret = c("specificity", "sensitivity"), transpose = FALSE))
  df <- data.frame(FPR = 1 - rc$specificity, TPR = rc$sensitivity); df <- df[order(df$FPR, df$TPR), ]
  ggplot(df, aes(FPR, TPR)) +
    geom_abline(linetype = 2, color = "grey60") + geom_path(linewidth = 1.1, color = "#117733") +
    coord_equal() +
    annotate("text", x = 0.97, y = 0.08, hjust = 1, label = sprintf("AUC = %.2f", auc_sel),
             color = "#117733", size = 4) +
    labs(x = "False positive rate", y = "True positive rate") + theme_fig
}
traj_plot <- function() {
  tr <- traj %>% mutate(Group = factor(Group, levels = c("Survivor", "NonSurvivor")),
                        Phase = factor(Phase, levels = c("Day18", "Week1", "Mid", "Late")))
  sm <- tr %>% group_by(Phase, Group) %>%
    summarise(mean = mean(Score), se = sd(Score) / sqrt(n()), .groups = "drop")
  ypos <- tr %>% group_by(Phase) %>% summarise(y = max(Score) + 0.05, .groups = "drop")
  ann <- left_join(trajp, ypos, by = "Phase") %>%
    mutate(Phase = factor(Phase, levels = c("Day18", "Week1", "Mid", "Late")))
  ggplot() +
    geom_line(data = sm, aes(Phase, mean, group = Group, color = Group), linewidth = 1) +
    geom_errorbar(data = sm, aes(Phase, ymin = mean - se, ymax = mean + se, color = Group), width = 0.18, linewidth = 0.5) +
    geom_point(data = tr, aes(Phase, Score, color = Group), position = position_jitter(width = 0.06, height = 0),
               alpha = 0.65, size = 1.5) +
    geom_point(data = sm, aes(Phase, mean, color = Group), size = 2.6) +
    geom_text(data = ann, aes(Phase, y, label = star), size = 5, inherit.aes = FALSE) +
    scale_color_manual(values = COLORS_BIN[c("Survivor", "NonSurvivor")]) +
    labs(x = NULL, y = "EN score (P[non-survival])") +
    theme_fig + theme(legend.position = "top", legend.title = element_blank())
}
null_plot <- function() {
  ggplot(data.frame(x = nulls$null_auc), aes(x)) +
    geom_histogram(bins = 25, fill = "grey78", color = "white") +
    geom_vline(xintercept = auc_sel, color = COLORS_BIN["NonSurvivor"], linewidth = 1.1) +
    annotate("text", x = auc_sel, y = Inf, vjust = 1.6, hjust = 1.05,
             label = sprintf("AUC %.2f\np= %.3f", auc_sel, perm_p), size = 3.2, color = COLORS_BIN["NonSurvivor"]) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(x = "Null AUC (permuted labels)", y = "Count") + theme_fig
}

## ---- house-style header ----
hdr <- function(letter, title) gTree(children = gList(
  textGrob(letter, x = unit(0.5, "lines"), y = 0.5, just = c("left", "center"), gp = gpar(fontsize = 16, fontface = "bold")),
  textGrob(title, x = 0.5, y = 0.5, just = "center", gp = gpar(fontsize = 12, fontface = "bold"))))
HH <- unit(1.5, "lines")

hmfile <- file.path(FIG_DIR, "_hm.png"); render_heatmap(hmfile)

panelA <- arrangeGrob(hdr("A", sprintf("Post-resection signature (%d genes)", n_genes)),
                      arrangeGrob(heatmap_grob(hmfile, "i.)"), pca_plot("ii.)"), ncol = 2),
                      heights = unit.c(HH, unit(1, "null")))
panelB <- arrangeGrob(hdr("B", "Cross-validated ROC"), roc_plot(), heights = unit.c(HH, unit(1, "null")))
panelC <- arrangeGrob(hdr("C", "EN-score trajectory"), traj_plot(), heights = unit.c(HH, unit(1, "null")))
panelE <- arrangeGrob(hdr("E", "Permutation validation"), null_plot(), heights = unit.c(HH, unit(1, "null")))
panelD_note <- arrangeGrob(hdr("D", "Scaffold ddPCR"),
                           textGrob("tdTomato copy number\n(see revision_figures/scripts/fig5_revised.R)",
                                    gp = gpar(fontsize = 10, col = "grey45", fontface = "italic")),
                           heights = unit.c(HH, unit(1, "null")))

## individual panels
ggsave(file.path(FIG_DIR, "5A_signature_pca.png"), pca_plot("ii.)"), width = 4.4, height = 3.4, dpi = 200)
ggsave(file.path(FIG_DIR, "5B_ROC.png"), roc_plot(), width = 3.6, height = 3.4, dpi = 200)
ggsave(file.path(FIG_DIR, "5C_trajectory.png"), traj_plot(), width = 5.0, height = 3.4, dpi = 200)
ggsave(file.path(FIG_DIR, "5E_permutation.png"), null_plot(), width = 4.2, height = 3.2, dpi = 200)

## composite (A row; B|C row; E|D-note row)
rowBC <- arrangeGrob(panelB, panelC, ncol = 2, widths = c(1, 1.25))
rowED <- arrangeGrob(panelE, panelD_note, ncol = 2, widths = c(1.4, 1))
fig <- arrangeGrob(panelA, rowBC, rowED, ncol = 1, heights = c(1.25, 1.0, 0.85))
png(file.path(FIG_DIR, "Figure5_postresection_housestyle.png"), width = 2500, height = 3100, res = 215)
grid.draw(fig); dev.off()
file.remove(hmfile)
message("Done: ", file.path(FIG_DIR, "Figure5_postresection_housestyle.png"))
