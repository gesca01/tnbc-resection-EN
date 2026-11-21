## Figure 5: A heatmap, B PCA, C leave-one-out cross-validated ROC, D EN-score trajectory,
## E tumor-derived tdTomato DNA (dPCR) at resection and week 1.
source("revision_figures/scripts/house_style.R")
source_upto(file.path(PROJ, "NCT_paper_pipeline/03_scripts_simplified/elastic_net/figure_postresection_signature_housestyle.R"))
suppressPackageStartupMessages(library(dplyr))
hs_apply()
HM_GENE_PT <- 7.5   # 30 genes
render_heatmap <- capture_heatmap(render_heatmap)


## Stars and 'ns' on one line above all points.
traj_plot <- function() {
  set.seed(1)
  lv <- c("Day18", "Week1", "Mid", "Late")
  tr <- traj %>% mutate(Group = factor(Group, levels = c("Survivor", "NonSurvivor")), Phase = factor(Phase, levels = lv))
  sm <- tr %>% group_by(Phase, Group) %>% summarise(mean = mean(Score), se = sd(Score) / sqrt(n()), .groups = "drop")
  ann <- trajp %>% mutate(Phase = factor(Phase, levels = lv), y = max(tr$Score) + 0.1)
  ggplot() +
    geom_line(data = sm, aes(Phase, mean, group = Group, color = Group), linewidth = 1) +
    geom_errorbar(data = sm, aes(Phase, ymin = mean - se, ymax = mean + se, color = Group), width = 0.18, linewidth = 0.5) +
    geom_point(data = tr, aes(Phase, Score, color = Group), position = position_jitter(width = 0.06, height = 0),
               alpha = 0.65, size = 1.5) +
    geom_point(data = sm, aes(Phase, mean, color = Group), size = 2.6) +
    geom_text(data = ann, aes(Phase, y, label = star), size = GEOM_PT(12), fontface = "bold", vjust = 0) +
    scale_y_continuous(limits = c(NA, max(tr$Score) + 0.2), breaks = seq(0, 1, 0.25)) +
    labs(x = NULL, y = "EN score (P[non-survival])")
}

pretty_phase <- c(Day18 = "Resection", Week1 = "Week 1", Mid = "Mid", Late = "Late")
pca_panel <- function() pca_plot(NULL) +
  scale_group_color() + scale_shape_manual(values = c(treated = 16, untreated = 17), labels = c("Treated", "Untreated")) +
  theme_house + theme(legend.position = "top", legend.title = element_blank(), legend.box = "vertical",
                      legend.spacing.y = unit(0, "pt"), legend.margin = margin(0, 0, 0, 0))
roc_panel <- function() roc_plot() + scale_x_continuous(breaks = seq(0, 1, 0.5)) +
  scale_y_continuous(breaks = seq(0, 1, 0.5)) + theme_house
traj_panel <- function() traj_plot() + scale_group_color() + scale_x_discrete(labels = pretty_phase) +
  labs(x = "Time after resection", y = "EN risk score\n(P[non-survival])") +
  theme_house + theme(legend.position = "top", legend.title = element_blank())

dpcr_timecourse <- function() {
  set.seed(1)
  d <- read_dpcr(); d <- d[d$timepoint != "Week 8", ]
  lv <- c("Resection", "Week 1"); d$x <- match(d$timepoint, lv)
  sm <- d %>% group_by(timepoint, group, x) %>%
    summarise(mean = mean(copies_per_uL), sd = sd(copies_per_uL), .groups = "drop")
  st <- data.frame(x = seq_along(lv), y = max(d$copies_per_uL) * 1.08,
                   lab = sapply(lv, function(t) p_stars(dpcr_test(d, t))))
  d$xj <- d$x + runif(nrow(d), -0.05, 0.05)
  ggplot() +
    geom_line(data = sm, aes(x, mean, color = group, group = group), linewidth = 0.9) +
    geom_errorbar(data = sm, aes(x, ymin = pmax(mean - sd, 0), ymax = mean + sd, color = group),
                  width = 0.12, linewidth = 0.5) +
    geom_point(data = d, aes(xj, copies_per_uL, color = group), size = 1.6, alpha = 0.7) +
    geom_point(data = sm, aes(x, mean, color = group), size = 2.8) +
    geom_text(data = st, aes(x, y, label = lab), size = GEOM_PT(14), fontface = "bold") +
    scale_group_color() +
    scale_x_continuous(breaks = seq_along(lv), labels = lv, expand = expansion(add = 0.4)) +
    labs(x = "Time after resection", y = DPCR_YLAB) +
    theme_house + theme(legend.position = "top", legend.title = element_blank())
}

render_heatmap("A")
pA <- wrap("A", sprintf("Post-resection signature (%d genes)", n_genes), hm_grob("A"))
pB <- wrap("B", "Principal component analysis", pca_panel())
pC <- wrap("C", "Cross-validated ROC", roc_panel())
pD <- wrap("D", "EN risk score over time", traj_panel())
pE <- wrap("E", "Tumor DNA in the niche (dPCR)", dpcr_timecourse())
layout5 <- function(a) arrangeGrob(a,
                   arrangeGrob(pB, pC, ncol = 2, widths = c(0.56, 0.44)),
                   arrangeGrob(pD, pE, ncol = 2, widths = c(0.56, 0.44)),
                   ncol = 1, heights = c(1.5, 1, 1))
save_fig(layout5(pA), "Fig5.png", height = 9.3)
d <- read_dpcr()
for (t in c("Resection", "Week 1")) message(sprintf("dPCR %s: Mann-Whitney p = %.4f", t, dpcr_test(d, t)))
