## Shared figure style for the manuscript figures (matched to Figure 1):
## Arial, bold axis titles and tick labels, black axes, bold panel letters and titles.
## Every figure is drawn at its final printed size (7 in wide, full page width),
## so the point sizes below are the sizes that appear in the manuscript.
source("revision_figures/scripts/common.R")
suppressPackageStartupMessages({ library(ggplot2); library(grid); library(gridExtra); library(pheatmap) })

FONT   <- "Arial"
PAGE_W <- 7.0                           # inches, full-width figure
PT <- list(letter = 18, title = 11, axis_title = 10, tick = 9, legend = 9, strip = 10,
           hm_gene = 8, hm_base = 9)
GEOM_PT <- function(pt) pt / .pt        # ggplot text size (mm) for a size in points
LINE_W  <- 0.6                          # axis / tick line width

theme_house <- theme_classic(base_size = PT$tick, base_family = FONT) +
  theme(axis.title   = element_text(size = PT$axis_title, face = "bold", color = "black"),
        axis.text    = element_text(size = PT$tick, face = "bold", color = "black"),
        axis.line    = element_line(linewidth = LINE_W, color = "black"),
        axis.ticks   = element_line(linewidth = LINE_W, color = "black"),
        axis.ticks.length = unit(3, "pt"),
        strip.background  = element_blank(),
        strip.text   = element_text(size = PT$strip, face = "bold", color = "black"),
        legend.text  = element_text(size = PT$legend, color = "black"),
        legend.title = element_text(size = PT$legend, face = "bold", color = "black"),
        legend.key.size = unit(0.9, "lines"),
        axis.title.y = element_text(margin = margin(r = 5)),
        plot.margin  = margin(4, 8, 4, 6))

## Panel header: bold letter at the left, bold centred title
hs_hdr <- function(letter, title = "") {
  lw <- unit(PT$letter + 8, "pt")   # letter column; title is centred in the space to its right
  gTree(children = gList(
    textGrob(letter, x = unit(2, "pt"), y = 0.5, just = c("left", "center"),
             gp = gpar(fontsize = PT$letter, fontface = "bold", fontfamily = FONT)),
    textGrob(title, x = lw + (unit(1, "npc") - lw) * 0.5, y = 0.5, just = "center",
             gp = gpar(fontsize = PT$title, fontface = "bold", fontfamily = FONT))))
}
HS_HH <- unit(PT$letter + 6, "pt")
wrap <- function(letter, title, g) arrangeGrob(hs_hdr(letter, title), g, heights = unit.c(HS_HH, unit(1, "null")))

## Call after source_upto() so the original panel builders pick up the house style
## Text is measured on the current device; keep an Arial PNG device open so layouts are
## computed with the real font metrics (the default PDF device has no Arial).
hs_measure_device <- function() if (names(dev.cur()) != "quartz_off_screen" && names(dev.cur()) != "png")
  png(file.path(tempdir(), "measure.png"), width = PAGE_W, height = 9, units = "in", res = 72, family = FONT)
hs_measure_device()
hs_apply <- function() {
  hs_measure_device()
  assign("theme_fig", theme_house, envir = globalenv())
  assign("hdr", hs_hdr, envir = globalenv()); assign("HH", HS_HH, envir = globalenv())
  assign("heatmap_grob", hm_grob, envir = globalenv())
}

## Heatmaps are kept as vector grobs (not pasted-in PNGs) so gene names print at their real size.
## The original scripts call render_heatmap(..., file) then heatmap_grob(file); wrapping the
## renderer with capture_heatmap() and using hm_grob() makes the file name just a key.
HM <- new.env()
pretty_levels <- c(NonSurvivor = "Non-survivor", Week1 = "Week 1")
pheatmap <- function(mat, ..., annotation_col = NULL, annotation_colors = NA,
                     fontsize_row = NULL, fontsize = NULL, silent = NULL) {
  if (!is.null(annotation_col)) for (k in names(annotation_col)) {
    f <- annotation_col[[k]]; lv <- levels(f); new <- ifelse(lv %in% names(pretty_levels), pretty_levels[lv], lv)
    annotation_col[[k]] <- factor(f, levels = lv, labels = new)
    if (is.list(annotation_colors) && !is.null(annotation_colors[[k]]))
      names(annotation_colors[[k]]) <- ifelse(names(annotation_colors[[k]]) %in% names(pretty_levels),
                                              pretty_levels[names(annotation_colors[[k]])], names(annotation_colors[[k]]))
  }
  ph <- pheatmap::pheatmap(mat, ..., annotation_col = annotation_col, annotation_colors = annotation_colors,
                           fontsize_row = HM_GENE_PT, fontsize = PT$hm_base, fontfamily = FONT,
                           treeheight_row = 20, treeheight_col = 15, silent = TRUE)
  assign(".last", ph$gtable, envir = HM); invisible(ph)
}
HM_GENE_PT <- PT$hm_gene
capture_heatmap <- function(render) { force(render); function(...) {
  a <- list(...); render(...); assign(a[[length(a)]], get(".last", envir = HM), envir = HM) }
}
hm_grob <- function(file, ...) get(file, envir = HM)

## Shared scales
lab_group <- c(Survivor = "Survivor", NonSurvivor = "Non-survivor")
scale_group_color <- function(...) scale_color_manual(values = c(Survivor = COL_S, NonSurvivor = COL_NS), labels = lab_group, ...)
scale_group_fill  <- function(...) scale_fill_manual(values = c(Survivor = COL_S, NonSurvivor = COL_NS), labels = lab_group, ...)

## Two-signature ROC (day 18 vs day 7) with a colour key that carries the AUCs
roc_two <- function(A, B) {
  rd <- function(d, lab) {
    rc <- as.data.frame(coords(d$roc, "all", ret = c("specificity", "sensitivity"), transpose = FALSE))
    out <- data.frame(FPR = 1 - rc$specificity, TPR = rc$sensitivity, Sig = lab); out[order(out$FPR, out$TPR), ]
  }
  labs_sig <- c(D18 = sprintf("Day 18, AUC = %.2f", A$auc), D7 = sprintf("Day 7, AUC = %.2f", B$auc))
  df <- rbind(rd(A, labs_sig[["D18"]]), rd(B, labs_sig[["D7"]]))
  df$Sig <- factor(df$Sig, levels = labs_sig)
  ggplot(df, aes(FPR, TPR, color = Sig)) +
    geom_abline(linetype = 2, color = "grey60") + geom_path(linewidth = 1.1) +
    scale_color_manual(values = setNames(c("#117733", "#882255"), labs_sig)) +
    scale_x_continuous(breaks = seq(0, 1, 0.5)) + scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    coord_equal() +
    labs(x = "False positive rate", y = "True positive rate", color = NULL) +
    theme_house + theme(legend.position = c(0.98, 0.04), legend.justification = c(1, 0),
                        legend.background = element_blank(), legend.key.width = unit(1.2, "lines"))
}

save_fig <- function(g, name, height, width = PAGE_W) {
  f <- file.path(OUT_DIR, name)
  cur <- dev.cur()
  png(f, width = width, height = height, units = "in", res = 300, family = FONT)
  grid.newpage(); grid.draw(g); dev.off(); dev.set(cur)
  message("wrote ", f)
}
