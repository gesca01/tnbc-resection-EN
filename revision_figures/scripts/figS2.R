## Figure S2: elastic-net signatures in untreated mice. Risk scores are the
## leave-one-out cross-validated scores used for the ROC.
source("revision_figures/scripts/house_style.R")
source_upto(file.path(PROJ, "Figure3_signature_analysis_report/scripts/19_suppfig_S1_noPCA_day18scores.R"))
hs_apply()
render_heatmap_clean <- capture_heatmap(render_heatmap_clean)

cv_scores_plot <- function() {
  sc <- rbind(data.frame(TP = "Resection day (day 18)", Score = A$scores, Group = A$group),
              data.frame(TP = "Pre-treatment (day 7)",  Score = B$scores, Group = B$group))
  sc$TP <- factor(sc$TP, levels = c("Resection day (day 18)", "Pre-treatment (day 7)"))
  sc$Group <- factor(sc$Group, levels = c("Survivor", "NonSurvivor"))
  pv <- sapply(split(sc, sc$TP), function(x) suppressWarnings(wilcox.test(Score ~ Group, data = x)$p.value))
  lab <- data.frame(TP = factor(names(pv), levels = levels(sc$TP)),
                    txt = ifelse(pv < 0.001, "Wilcoxon p < 0.001", sprintf("Wilcoxon p = %.3f", pv)))
  ggplot(sc, aes(Group, Score, fill = Group)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.6) +
    geom_jitter(aes(color = Group), width = 0.12, size = 2, alpha = 0.9) +
    geom_text(data = lab, aes(x = 1.5, y = 1.08, label = txt), inherit.aes = FALSE, size = GEOM_PT(9)) +
    facet_wrap(~TP) +
    scale_fill_manual(values = COLORS_BIN) + scale_color_manual(values = COLORS_BIN) +
    scale_x_discrete(labels = c(Survivor = "Survivor", NonSurvivor = "Non-survivor")) +
    coord_cartesian(ylim = c(0, 1.12)) +
    labs(x = NULL, y = "Cross-validated EN score\n(predicted probability of non-survival)") +
    theme_fig + theme(legend.position = "none")
}
## STYLE=dots: same cross-validated scores as dots + median bar + 0.5 decision line (matches the Fig 3E dots option)
STYLE <- Sys.getenv("STYLE", "box")
cv_scores_dots <- function() {
  sc <- rbind(data.frame(TP = "Resection day (day 18)", Score = A$scores, Group = A$group),
              data.frame(TP = "Pre-treatment (day 7)",  Score = B$scores, Group = B$group))
  sc$TP <- factor(sc$TP, levels = c("Resection day (day 18)", "Pre-treatment (day 7)"))
  sc$Group <- factor(sc$Group, levels = c("Survivor", "NonSurvivor"))
  pv <- sapply(split(sc, sc$TP), function(x) suppressWarnings(wilcox.test(Score ~ Group, data = x)$p.value))
  lab <- data.frame(TP = factor(names(pv), levels = levels(sc$TP)),
                    txt = ifelse(pv < 0.001, "p < 0.001", sprintf("p = %.3f", pv)))
  ggplot(sc, aes(Group, Score, color = Group)) +
    geom_hline(yintercept = 0.5, linetype = 2, color = "grey60") +
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.45, linewidth = 0.4, color = "grey20") +
    geom_jitter(width = 0.08, height = 0, size = 2.8, alpha = 0.9) +
    geom_text(data = lab, aes(x = 1.5, y = 1.1, label = txt), inherit.aes = FALSE, size = GEOM_PT(9)) +
    facet_wrap(~TP) +
    scale_color_manual(values = COLORS_BIN) +
    scale_x_discrete(labels = c(Survivor = "Survivor", NonSurvivor = "Non-survivor")) +
    scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25)) +
    labs(x = NULL, y = "Cross-validated EN risk score\n(predicted probability of non-survival)") +
    theme_fig + theme(legend.position = "none")
}
perm_plot <- function(d) null_plot_clean(d, NULL, NULL) +
  labs(x = "Cross-validated AUC with shuffled labels", y = "Number of permutations")

render_heatmap_clean(A, "A"); render_heatmap_clean(B, "B")
pA <- wrap("A", sprintf("Resection-day signature (%d genes)", length(A$sig_genes)), hm_grob("A"))
pB <- wrap("B", sprintf("Pre-treatment signature (%d genes)", length(B$sig_genes)), hm_grob("B"))
pC <- wrap("C", "Cross-validated ROC", roc_two(A, B))
pD <- wrap("D", "Cross-validated EN risk scores", (if (STYLE == "dots") cv_scores_dots() else cv_scores_plot()) +
             scale_y_continuous(breaks = seq(0, 1, 0.25)) + labs(y = "Cross-validated EN risk score\n(P[non-survival])"))
pE <- wrap("E", "Permutation test, resection day", perm_plot(A) + labs(x = "AUC with shuffled labels", y = "Permutations"))
pF <- wrap("F", "Permutation test, pre-treatment", perm_plot(B) + labs(x = "AUC with shuffled labels", y = "Permutations"))
fig <- arrangeGrob(arrangeGrob(pA, pB, ncol = 2), arrangeGrob(pC, pD, ncol = 2, widths = c(0.42, 0.58)),
                   arrangeGrob(pE, pF, ncol = 2), ncol = 1, heights = c(1.2, 1.05, 0.8))
save_fig(fig, if (STYLE == "dots") "FigS2_dots.png" else "FigS2.png", height = 8.8)
cat(sprintf("Day 18: n=%d (%d S / %d NS), %d genes, AUC %.2f (%.2f-%.2f), perm p=%.3f\n", A$n,
            sum(A$group == "Survivor"), sum(A$group == "NonSurvivor"), length(A$sig_genes), A$auc, A$ci[1], A$ci[3], A$perm_p))
cat(sprintf("Day 7: n=%d (%d S / %d NS), %d genes, AUC %.2f (%.2f-%.2f), perm p=%.3f\n", B$n,
            sum(B$group == "Survivor"), sum(B$group == "NonSurvivor"), length(B$sig_genes), B$auc, B$ci[1], B$ci[3], B$perm_p))
