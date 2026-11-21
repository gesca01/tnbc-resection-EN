## Figure 3: A dPCR at resection, B day-18 heatmap, C day-7 heatmap,
## D leave-one-out cross-validated ROC, E elastic-net risk scores.
## Writes two versions that differ only in panel E:
##   Fig3.png                             scores from the model fit on all treated mice (manuscript)
##   Fig3_alt_crossvalidated_scores.png   leave-one-out scores (the same scores the ROC in D is built from)
source("revision_figures/scripts/house_style.R")
SIG_ROOT <- file.path(PROJ, "NCT_paper_pipeline/02_cache/elastic_net_day7_day18_source_match/results_EN_signatures")
sig_all <- unlist(lapply(c("Day18_Treated_SvNS", "Day7_Treated_SvNS"),
                         function(k) read.csv(file.path(SIG_ROOT, k, "selected_signature_genes.csv"))$gene))
suppressPackageStartupMessages(library(dplyr))
source_upto(file.path(PROJ, "Figure3_signature_analysis_report/scripts/10_figure3_housestyle.R"),
            keep_genes = sig_all)
stopifnot(length(A$sig) == 21, length(B$sig) == 15)
hs_apply()
render_heatmap <- capture_heatmap(render_heatmap)
refit_scores <- list(A = A$refit, B = B$refit)
cv_scores    <- list(A = A$roc$predictor, B = B$roc$predictor)
stopifnot(length(cv_scores$A) == A$n, length(cv_scores$B) == B$n)

dpcr_bar <- function() {
  d <- read_dpcr(); d <- d[d$timepoint == "Resection", ]
  sm <- d %>% group_by(group) %>% summarise(m = mean(copies_per_uL), s = sd(copies_per_uL), .groups = "drop")
  ytop <- max(d$copies_per_uL) * 1.1
  ggplot() +
    geom_col(data = sm, aes(group, m, fill = group), width = 0.6, alpha = 0.45, color = "black", linewidth = LINE_W) +
    geom_errorbar(data = sm, aes(group, ymin = pmax(m - s, 0), ymax = m + s), width = 0.2, linewidth = LINE_W) +
    geom_jitter(data = d, aes(group, copies_per_uL, color = group), width = 0.1, height = 0, size = 2) +
    annotate("segment", x = 1, xend = 2, y = ytop, yend = ytop, linewidth = LINE_W) +
    annotate("text", x = 1.5, y = ytop * 1.02, vjust = 0, size = GEOM_PT(14), fontface = "bold",
             label = p_stars(dpcr_test(read_dpcr(), "Resection"))) +
    scale_group_fill() + scale_group_color() +
    scale_x_discrete(labels = c(Survivor = "Survivors", NonSurvivor = "Non-survivors")) +
    scale_y_continuous(limits = c(0, ytop * 1.12), expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = DPCR_YLAB) + theme_house + theme(legend.position = "none")
}

scores_panel <- function(sc_list, cv) {
  sc <- rbind(data.frame(TP = "Day 7 (pre-treatment)", Score = sc_list$B, Group = B$group),
              data.frame(TP = "Day 18 (resection)",    Score = sc_list$A, Group = A$group))
  sc$TP <- factor(sc$TP, levels = c("Day 7 (pre-treatment)", "Day 18 (resection)"))
  sc$Group <- factor(sc$Group, levels = c("Survivor", "NonSurvivor"))
  pv <- sapply(split(sc, sc$TP), function(x) wilcox.test(Score ~ Group, data = x)$p.value)
  nn <- sapply(split(sc, sc$TP), function(x) sprintf("n = %d vs %d", sum(x$Group == "Survivor"), sum(x$Group != "Survivor")))
  lab <- data.frame(TP = factor(names(pv), levels = levels(sc$TP)),
                    txt = paste0(ifelse(pv < 0.001, "p < 0.001", sprintf("p = %.3f", pv)), "\n", nn))
  ggplot(sc, aes(Group, Score)) +
    geom_boxplot(aes(fill = Group), width = 0.55, outlier.shape = NA, alpha = 0.45, linewidth = 0.5) +
    geom_jitter(aes(color = Group), width = 0.1, height = 0, size = 1.8) +
    geom_text(data = lab, aes(x = 1.5, y = if (cv) 1.13 else 0.5, label = txt), size = GEOM_PT(9), lineheight = 0.9) +
    facet_wrap(~TP) +
    scale_group_fill() + scale_group_color() +
    scale_x_discrete(labels = c(Survivor = "Survivor", NonSurvivor = "Non-survivor")) +
    scale_y_continuous(limits = c(0, if (cv) 1.22 else 1.02), breaks = seq(0, 1, 0.25)) +
    labs(x = NULL, y = if (cv) "Cross-validated EN risk score\n(P[non-survival])" else "EN risk score\n(P[non-survival])") +
    theme_house + theme(legend.position = "none")
}

render_heatmap(A, "A"); render_heatmap(B, "B")
pA <- wrap("A", "Tumor DNA at resection", dpcr_bar())
pB <- wrap("B", sprintf("Resection-day signature (%d genes)", length(A$sig)), hm_grob("A"))
pC <- wrap("C", sprintf("Pre-treatment signature (%d genes)", length(B$sig)), hm_grob("B"))
pD <- wrap("D", "Cross-validated ROC", roc_two(A, B))
build <- function(cv, name) {
  pE <- wrap("E", if (cv) "Cross-validated EN risk scores" else "EN risk scores (model fit on all mice)",
             scores_panel(if (cv) cv_scores else refit_scores, cv))
  fig <- arrangeGrob(arrangeGrob(pA, pB, ncol = 2, widths = c(0.38, 0.62)),
                     arrangeGrob(pC, pD, ncol = 2, widths = c(0.58, 0.42)),
                     arrangeGrob(pE, nullGrob(), ncol = 2, widths = c(0.62, 0.38)),
                     ncol = 1, heights = c(1.15, 1.05, 0.85))
  save_fig(fig, name, height = 9)
}
build(FALSE, "Fig3.png")
build(TRUE,  "Fig3_alt_crossvalidated_scores.png")
for (k in c("A", "B")) {
  d <- get(k); s <- list(refit = refit_scores[[k]], cv = cv_scores[[k]])
  for (t in names(s)) message(sprintf("%s %s: Wilcoxon p = %.4f", if (k == "A") "Day18" else "Day7", t,
                                      wilcox.test(s[[t]][d$group == "Survivor"], s[[t]][d$group != "Survivor"])$p.value))
}
cat(sprintf("D18 AUC %.3f | D7 AUC %.3f\n", A$auc, B$auc))
