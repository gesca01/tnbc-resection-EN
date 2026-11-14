## Figure 2.
## A: tumor volume change, day 7 to day 17 (data/fig2A_tumor_percent_change.csv)
## B: tumor volume change vs survival time, treated mice (data/fig2B_survival_vs_tumor_change.csv)
## C: representative IVIS image (data/fig2C_ivis_image.png)
## D: lung IVIS signal one day after resection (counts)

source("revision_figures/scripts/common.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(grid); library(gridExtra) })
set.seed(1)  # fixed jitter so the figure is reproducible
theme_fig <- theme_classic(base_size = 12) + theme(plot.margin = margin(6, 8, 6, 8))
COLS <- c(Survivor = COL_S, NonSurvivor = COL_NS)
GLAB <- c(Survivor = "Survivors", NonSurvivor = "Non-survivors")

d <- read.csv(file.path(WD, "data", "fig2A_tumor_percent_change.csv"))
d$group <- factor(d$group, levels = c("Survivor", "NonSurvivor"))
d$cohort <- factor(d$cohort, levels = c("Treated", "Untreated"))
sm <- d %>% group_by(cohort, group) %>% summarise(m = mean(percent_change_D7_D17), s = sd(percent_change_D7_D17), n = n(), .groups = "drop")
pv <- sapply(split(d, d$cohort), function(x) suppressWarnings(wilcox.test(percent_change_D7_D17 ~ group, data = x)$p.value))
lab <- data.frame(cohort = factor(names(pv), levels = levels(d$cohort)),
                  txt = ifelse(pv < 0.05, sprintf("p = %.3f", pv), sprintf("ns (p = %.2f)", pv)))
pA <- ggplot() +
  geom_col(data = sm, aes(group, m, fill = group), width = 0.6, alpha = 0.5, color = "black") +
  geom_errorbar(data = sm, aes(group, ymin = m - s, ymax = m + s), width = 0.2) +
  geom_jitter(data = d, aes(group, percent_change_D7_D17, color = group), width = 0.1, height = 0, size = 2.2) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(data = lab, aes(x = 1.5, y = 1000, label = txt), size = 3.6) +
  facet_wrap(~cohort) +
  scale_fill_manual(values = COLS) + scale_color_manual(values = COLS) +
  scale_x_discrete(labels = GLAB) +
  scale_y_continuous(limits = c(-100, 1050), breaks = seq(0, 1000, 250)) +
  labs(x = NULL, y = "Primary tumor volume change (%)\n(Day 7 to Day 17)") +
  theme_fig + theme(legend.position = "none")

ivis <- data.frame(group = factor(rep(c("Survivor", "NonSurvivor"), each = 2), levels = c("Survivor", "NonSurvivor")),
                   signal = factor(rep(c("Lung signal", "No signal"), 2), levels = c("Lung signal", "No signal")),
                   n = c(3, 4, 4, 4))
pD <- ggplot(ivis, aes(group, n, fill = signal)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, color = "black") +
  geom_text(aes(label = n), position = position_dodge(width = 0.75), vjust = -0.4, size = 3.6) +
  scale_fill_manual(values = c("Lung signal" = "grey25", "No signal" = "grey85")) +
  scale_x_discrete(labels = GLAB) +
  scale_y_continuous(limits = c(0, 5.2), breaks = 0:5, expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Number of mice", fill = "IVIS one day\nafter resection") +
  theme_fig + theme(legend.position = "right")

b <- read.csv(file.path(WD, "data", "fig2B_survival_vs_tumor_change.csv"))
ns <- subset(b, group == "NonSurvivor"); sv <- subset(b, group == "Survivor")
fit <- lm(percent_change_D7_D17 ~ survival_weeks, data = ns)
r2 <- summary(fit)$r.squared
SURV_X <- 8.3
pB <- ggplot() +
  annotate("rect", xmin = 7.55, xmax = 9.05, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_smooth(data = ns, aes(survival_weeks, percent_change_D7_D17), method = "lm", formula = y ~ x,
              se = FALSE, color = COL_NS, linewidth = 0.6, linetype = "dashed") +
  geom_point(data = ns, aes(survival_weeks, percent_change_D7_D17), color = COL_NS, shape = 16, size = 2.4) +
  geom_jitter(data = sv, aes(SURV_X, percent_change_D7_D17), color = COL_S, shape = 17, size = 2.4,
              width = 0.18, height = 0) +
  annotate("text", x = 2, y = 380, label = sprintf("Non-survivors: R\u00b2 = %.3f", r2),
           hjust = 0, size = 3.6) +
  scale_x_continuous(limits = c(2, 9.05), breaks = c(2, 3, 4, 5, 6, 7, SURV_X),
                     labels = c(2:7, "Survivors"), expand = expansion(mult = c(0.02, 0))) +
  scale_y_continuous(limits = c(-100, 400), breaks = seq(-100, 400, 100)) +
  labs(x = "Survival time after resection (weeks)",
       y = "Primary tumor volume change (%)\n(Day 7 to Day 17)") +
  theme_fig + theme(axis.text.x = element_text(face = c(rep("plain", 6), "bold")))

orig <- png::readPNG(file.path(WD, "data", "fig2C_ivis_image.png"))
pC <- rasterGrob(orig[205:648, 800:1420, 1:3], interpolate = TRUE)

hdr <- function(letter) textGrob(letter, x = unit(0.3, "lines"), just = "left", gp = gpar(fontsize = 16, fontface = "bold"))
for (k in list(list("A", pA, "Fig2A_rebuilt.png", 6, 4.2), list("D", pD, "Fig2D_rebuilt.png", 5, 3.6))) {
  png(file.path(OUT_DIR, k[[3]]), width = k[[4]], height = k[[5]], units = "in", res = 300)
  grid.draw(arrangeGrob(hdr(k[[1]]), k[[2]], heights = unit.c(unit(1.6, "lines"), unit(1, "null")))); dev.off()
}
for (k in list(list("B", pB, "Fig2B_rebuilt.png", 6, 4.2))) {
  png(file.path(OUT_DIR, k[[3]]), width = k[[4]], height = k[[5]], units = "in", res = 300)
  grid.draw(arrangeGrob(hdr(k[[1]]), k[[2]], heights = unit.c(unit(1.6, "lines"), unit(1, "null")))); dev.off()
}
panel <- function(letter, g) arrangeGrob(hdr(letter), g, heights = unit.c(unit(1.6, "lines"), unit(1, "null")))
png(file.path(OUT_DIR, "Fig2_rebuilt.png"), width = 12, height = 8.8, units = "in", res = 300)
grid.draw(arrangeGrob(panel("A", pA), panel("C", pC), panel("B", pB), panel("D", pD),
                      ncol = 2, widths = c(1.05, 1), layout_matrix = rbind(c(1, 2), c(3, 4))))
dev.off()
print(sm); print(pv); cat("Fig 2B R2 =", r2, "\n")
