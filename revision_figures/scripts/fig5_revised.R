## Figure 5: A heatmap, B PCA, C ROC, D EN-score trajectory,
## E tumor-derived tdTomato DNA (dPCR) at resection and week 1.

source("revision_figures/scripts/common.R")
source_upto(file.path(PROJ, "NCT_paper_pipeline/03_scripts_simplified/elastic_net/figure_postresection_signature_housestyle.R"))
library(dplyr)

pretty_phase <- c(Day18 = "Resection", Week1 = "Week 1", Mid = "Mid", Late = "Late")
traj_plot_labeled <- function() traj_plot() + scale_x_discrete(labels = pretty_phase)

dpcr_timecourse <- function(include_week8) {
  d <- read_dpcr()
  if (!include_week8) d <- d[d$timepoint != "Week 8", ]
  lv <- unique(d$timepoint); d$x <- match(d$timepoint, lv)
  sm <- d %>% group_by(timepoint, group, x) %>%
    summarise(mean = mean(copies_per_uL), sd = sd(copies_per_uL), .groups = "drop")
  tps <- intersect(c("Resection", "Week 1"), lv)
  st <- data.frame(x = match(tps, lv), y = max(d$copies_per_uL) * 1.08,
                   lab = sapply(tps, function(t) p_stars(dpcr_test(d, t))))
  dodge <- c(Survivor = 0, NonSurvivor = 0)  # both groups on the timepoint tick
  d$xj <- d$x + dodge[as.character(d$group)] + runif(nrow(d), -0.05, 0.05)
  sm$xd <- sm$x + dodge[as.character(sm$group)]
  ggplot() +
    geom_line(data = sm, aes(xd, mean, color = group, group = group), linewidth = 0.9) +
    geom_errorbar(data = sm, aes(xd, ymin = pmax(mean - sd, 0), ymax = mean + sd, color = group),
                  width = 0.15, linewidth = 0.5) +
    geom_point(data = d, aes(xj, copies_per_uL, color = group, shape = group), size = 2, alpha = 0.8) +
    geom_point(data = sm, aes(xd, mean, color = group), size = 3) +
    geom_text(data = st, aes(x, y, label = lab), size = 5) +
    scale_color_manual(values = c(Survivor = COL_S, NonSurvivor = COL_NS), labels = c("Survivor", "Non-survivor")) +
    scale_shape_manual(values = c(Survivor = 16, NonSurvivor = 17), labels = c("Survivor", "Non-survivor")) +
    scale_x_continuous(breaks = seq_along(lv), labels = lv, expand = expansion(add = 0.4)) +
    labs(x = "Time after resection", y = DPCR_YLAB) +
    theme_fig + theme(legend.position = "top", legend.title = element_blank())
}

build <- function(include_week8, file) {
  hmfile <- file.path(TMP_DIR, "_hm5.png"); render_heatmap(hmfile)
  pA <- arrangeGrob(hdr("A", sprintf("Post-resection signature (%d genes)", n_genes)),
                    heatmap_grob(hmfile, ""), heights = unit.c(HH, unit(1, "null")))
  pB <- arrangeGrob(hdr("B", "Principal component analysis"), pca_plot(NULL), heights = unit.c(HH, unit(1, "null")))
  pC <- arrangeGrob(hdr("C", "Cross-validated ROC"), roc_plot(), heights = unit.c(HH, unit(1, "null")))
  pD <- arrangeGrob(hdr("D", "EN-score trajectory"), traj_plot_labeled(), heights = unit.c(HH, unit(1, "null")))
  pE <- arrangeGrob(hdr("E", "Scaffold tumor-derived DNA (digital PCR)"), dpcr_timecourse(include_week8),
                    heights = unit.c(HH, unit(1, "null")))
  fig <- arrangeGrob(arrangeGrob(pA, pB, ncol = 2, widths = c(1.15, 1)),
                     arrangeGrob(pC, pD, ncol = 2, widths = c(1, 1.25)),
                     arrangeGrob(pE, nullGrob(), ncol = 2, widths = c(1.2, 1)),
                     ncol = 1, heights = c(1.2, 1, 0.95))
  png(file, width = 2500, height = 3100, res = 215); grid.draw(fig); dev.off()
  message("wrote ", file)
}
build(FALSE, file.path(OUT_DIR, "Fig5_revised.png"))
build(TRUE,  file.path(OUT_DIR, "Fig5_revised_with_week8.png"))
d <- read_dpcr()
for (t in c("Resection", "Week 1"))
  message(sprintf("dPCR %s: Mann-Whitney p = %.4f", t, dpcr_test(d, t)))
