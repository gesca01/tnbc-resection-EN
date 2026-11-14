## Figure 3: A dPCR at resection, B day-18 heatmap, C day-7 heatmap,
## D cross-validated ROC, E elastic-net risk scores.

source("revision_figures/scripts/common.R")
SIG_ROOT <- file.path(PROJ, "NCT_paper_pipeline/02_cache/elastic_net_day7_day18_source_match/results_EN_signatures")
sig_all <- unlist(lapply(c("Day18_Treated_SvNS", "Day7_Treated_SvNS"),
                         function(k) read.csv(file.path(SIG_ROOT, k, "selected_signature_genes.csv"))$gene))
suppressPackageStartupMessages(library(dplyr))
source_upto(file.path(PROJ, "Figure3_signature_analysis_report/scripts/10_figure3_housestyle.R"),
            keep_genes = sig_all)
stopifnot(length(A$sig) == 21, length(B$sig) == 15)

dpcr_bar <- function() {
  d <- read_dpcr(); d <- d[d$timepoint == "Resection", ]
  sm <- aggregate(copies_per_uL ~ group, d, function(v) c(m = mean(v), s = sd(v)))
  sm <- data.frame(group = sm$group, m = sm$copies_per_uL[, "m"], s = sm$copies_per_uL[, "s"])
  ytop <- max(d$copies_per_uL) * 1.12
  ggplot() +
    geom_col(data = sm, aes(group, m, fill = group), width = 0.6, alpha = 0.55, color = "black") +
    geom_errorbar(data = sm, aes(group, ymin = pmax(m - s, 0), ymax = m + s), width = 0.2) +
    geom_jitter(data = d, aes(group, copies_per_uL, color = group), width = 0.1, height = 0, size = 2.2) +
    annotate("segment", x = 1, xend = 2, y = ytop, yend = ytop) +
    annotate("text", x = 1.5, y = ytop * 1.03, label = p_stars(dpcr_test(read_dpcr(), "Resection")), size = 6) +
    scale_fill_manual(values = COLORS_BIN) + scale_color_manual(values = COLORS_BIN) +
    scale_x_discrete(labels = c(Survivor = "Survivors", NonSurvivor = "Non-survivors")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = DPCR_YLAB) + theme_fig + theme(legend.position = "none")
}
## SCORES=cv renders panel E from the leave-one-out scores (the ROC's own predictions) instead of refit-on-all
## SCORES=cv_dot / cv_waterfall: same cross-validated scores, cleaner layouts
STYLE  <- Sys.getenv("SCORES", "refit")
USE_CV <- STYLE %in% c("cv", "cv_dot", "cv_waterfall")
if (USE_CV) { A$refit <- A$roc$predictor; B$refit <- B$roc$predictor }
scores_with_p <- function() {
  pv <- c("Day 7 (pre-treatment)"  = wilcox.test(B$refit[B$group == "Survivor"], B$refit[B$group != "Survivor"])$p.value,
          "Day 18 (resection)"     = wilcox.test(A$refit[A$group == "Survivor"], A$refit[A$group != "Survivor"])$p.value)
  nn <- c("Day 7 (pre-treatment)" = sprintf("n = %d vs %d", sum(B$group == "Survivor"), sum(B$group != "Survivor")),
          "Day 18 (resection)"    = sprintf("n = %d vs %d", sum(A$group == "Survivor"), sum(A$group != "Survivor")))
  lab <- data.frame(TP = factor(names(pv), levels = c("Day 7 (pre-treatment)", "Day 18 (resection)")),
                    txt = paste0(ifelse(pv < 0.001, "p < 0.001", sprintf("p = %.3f", pv)), "\n", nn[names(pv)]))
  scores_plot() + geom_text(data = lab, aes(x = 1.5, y = if (USE_CV) 1.1 else 0.5, label = txt), inherit.aes = FALSE, size = 3.3) +
    scale_x_discrete(labels = c(Survivor = "Survivor", NonSurvivor = "Non-survivor")) +
    labs(y = if (USE_CV) "Cross-validated EN risk score\n(predicted probability of non-survival)"
             else "EN risk score\n(predicted probability of non-survival)")
}

tmpf <- function(n) file.path(TMP_DIR, paste0("_hm3_", n, ".png"))
cv_df <- function() {
  sc <- rbind(data.frame(TP = "Day 7 (pre-treatment)", Score = B$refit, Group = B$group),
              data.frame(TP = "Day 18 (resection)", Score = A$refit, Group = A$group))
  sc$TP <- factor(sc$TP, levels = c("Day 7 (pre-treatment)", "Day 18 (resection)"))
  sc$Group <- factor(sc$Group, levels = c("Survivor", "NonSurvivor")); sc
}
cv_p_labels <- function(sc) {
  pv <- sapply(split(sc, sc$TP), function(x) wilcox.test(Score ~ Group, data = x)$p.value)
  data.frame(TP = factor(names(pv), levels = levels(sc$TP)),
             txt = ifelse(pv < 0.001, "p < 0.001", sprintf("p = %.3f", pv)))
}
## Option 1: dots with median bars and the 0.5 decision line (no boxes)
scores_dot <- function() {
  sc <- cv_df(); lab <- cv_p_labels(sc)
  ggplot(sc, aes(Group, Score, color = Group)) +
    geom_hline(yintercept = 0.5, linetype = 2, color = "grey60") +
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.45, linewidth = 0.4, color = "grey20") +
    geom_jitter(width = 0.08, height = 0, size = 2.8, alpha = 0.9) +
    geom_text(data = lab, aes(x = 1.5, y = 1.1, label = txt), inherit.aes = FALSE, size = 3.4) +
    facet_wrap(~TP) +
    scale_color_manual(values = COLORS_BIN) +
    scale_x_discrete(labels = c(Survivor = "Survivor", NonSurvivor = "Non-survivor")) +
    scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25)) +
    labs(x = NULL, y = "Cross-validated EN risk score\n(predicted probability of non-survival)") +
    theme_fig + theme(legend.position = "none")
}
## Option 2: waterfall, each bar is one mouse ranked by score, colored by its actual outcome
scores_waterfall <- function() {
  sc <- cv_df() %>% group_by(TP) %>% arrange(Score, .by_group = TRUE) %>% mutate(rank = row_number()) %>% ungroup()
  lab <- cv_p_labels(sc)
  ggplot(sc, aes(rank, Score, fill = Group)) +
    geom_col(width = 0.8) +
    geom_hline(yintercept = 0.5, linetype = 2, color = "grey40") +
    geom_text(data = lab, aes(x = 1, y = 1.08, label = txt), inherit.aes = FALSE, hjust = 0, size = 3.3) +
    facet_wrap(~TP, scales = "free_x") +
    scale_fill_manual(values = COLORS_BIN, labels = c("Survivor", "Non-survivor")) +
    scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0))) +
    labs(x = "Mice ranked by score", y = "Cross-validated EN risk score\n(predicted probability of non-survival)") +
    theme_fig + theme(legend.position = "top", legend.title = element_blank(),
                      axis.text.x = element_blank(), axis.ticks.x = element_blank())
}
render_heatmap(A, tmpf("A")); render_heatmap(B, tmpf("B"))
wrap <- function(letter, title, g) arrangeGrob(hdr(letter, title), g, heights = unit.c(HH, unit(1, "null")))
pA <- wrap("A", "Scaffold tumor-derived DNA at resection", dpcr_bar())
pB <- wrap("B", sprintf("Resection-day signature (%d genes)", length(A$sig)), heatmap_grob(tmpf("A"), ""))
pC <- wrap("C", sprintf("Pre-treatment signature (%d genes)", length(B$sig)), heatmap_grob(tmpf("B"), ""))
pD <- wrap("D", "Cross-validated ROC", roc_plot())
pE <- wrap("E", if (USE_CV) "Cross-validated EN risk scores" else "EN risk scores",
           switch(STYLE, cv_dot = scores_dot(), cv_waterfall = scores_waterfall(), scores_with_p()))
fig <- arrangeGrob(arrangeGrob(pA, pB, ncol = 2, widths = c(0.75, 1.25)),
                   arrangeGrob(pC, pD, ncol = 2, widths = c(1.25, 0.9)),
                   arrangeGrob(pE, nullGrob(), ncol = 2, widths = c(1.2, 0.8)),
                   ncol = 1, heights = c(1.1, 1.1, 0.9))
out <- file.path(OUT_DIR, switch(STYLE, cv = "Fig3_revised_CVscores.png", cv_dot = "Fig3_revised_CVscores_dots.png",
                                 cv_waterfall = "Fig3_revised_CVscores_waterfall.png", "Fig3_revised.png"))
png(out, width = 2500, height = 3000, res = 215); grid.draw(fig); dev.off()
cat(sprintf("D18 AUC %.3f perm p %.3f | D7 AUC %.3f perm p %.3f\n", A$auc, A$perm_p, B$auc, B$perm_p))
