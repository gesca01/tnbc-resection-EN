## Figure 2 (drawn in the manuscript house style, see house_style.R).
## A: tumor volume change, day 7 to day 17 (data/fig2A_tumor_percent_change.csv)
## B: tumor volume change vs survival time, treated mice (data/fig2B_survival_vs_tumor_change.csv)
## C: representative IVIS image (data/fig2C_ivis_image.png)
## D: lung IVIS signal one day after resection (counts)
source("revision_figures/scripts/house_style.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(grid); library(gridExtra) })
set.seed(1)  # fixed jitter so the figure is reproducible
theme_fig <- theme_house
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
  geom_col(data = sm, aes(group, m, fill = group), width = 0.6, alpha = 0.5, color = "black", linewidth = LINE_W) +
  geom_errorbar(data = sm, aes(group, ymin = m - s, ymax = m + s), width = 0.2, linewidth = LINE_W) +
  geom_jitter(data = d, aes(group, percent_change_D7_D17, color = group), width = 0.1, height = 0, size = 1.8) +
  geom_hline(yintercept = 0, linewidth = LINE_W) +
  geom_text(data = lab, aes(x = 1.5, y = 1000, label = txt), size = GEOM_PT(9)) +
  facet_wrap(~cohort) +
  scale_fill_manual(values = COLS) + scale_color_manual(values = COLS) +
  scale_x_discrete(labels = c(Survivor = "Survivors", NonSurvivor = "Non-\nsurvivors")) +
  scale_y_continuous(limits = c(-100, 1050), breaks = seq(0, 1000, 250)) +
  labs(x = NULL, y = "Primary tumor volume change (%)\n(Day 7 to Day 17)") +
  theme_fig + theme(legend.position = "none")

ivis <- data.frame(group = factor(rep(c("Survivor", "NonSurvivor"), each = 2), levels = c("Survivor", "NonSurvivor")),
                   signal = factor(rep(c("Lung signal", "No signal"), 2), levels = c("Lung signal", "No signal")),
                   n = c(3, 4, 4, 4))
pD <- ggplot(ivis, aes(group, n, fill = signal)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, color = "black") +
  geom_text(aes(label = n), position = position_dodge(width = 0.75), vjust = -0.4, size = GEOM_PT(9)) +
  scale_fill_manual(values = c("Lung signal" = "grey25", "No signal" = "grey85")) +
  scale_x_discrete(labels = GLAB) +
  scale_y_continuous(limits = c(0, 6.2), breaks = 0:6, expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Number of mice", fill = NULL) +
  theme_fig + theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1))

b <- read.csv(file.path(WD, "data", "fig2B_survival_vs_tumor_change.csv"))
ns <- subset(b, group == "NonSurvivor"); sv <- subset(b, group == "Survivor")
fit <- lm(percent_change_D7_D17 ~ survival_weeks, data = ns)
r2 <- summary(fit)$r.squared
SURV_X <- 8.3
pB <- ggplot() +
  annotate("rect", xmin = 7.55, xmax = 9.05, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_hline(yintercept = 0, linewidth = LINE_W) +
  geom_smooth(data = ns, aes(survival_weeks, percent_change_D7_D17), method = "lm", formula = y ~ x,
              se = FALSE, color = COL_NS, linewidth = 0.6, linetype = "dashed") +
  geom_point(data = ns, aes(survival_weeks, percent_change_D7_D17), color = COL_NS, shape = 16, size = 2) +
  geom_jitter(data = sv, aes(SURV_X, percent_change_D7_D17), color = COL_S, shape = 17, size = 2.4,
              width = 0.18, height = 0) +
  annotate("text", x = 2, y = 380, label = sprintf("Non-survivors: R\u00b2 = %.3f", r2),
           hjust = 0, size = GEOM_PT(9)) +
  scale_x_continuous(limits = c(2, 9.05), breaks = c(2, 3, 4, 5, 6, 7, SURV_X),
                     labels = c(2:7, "Survivors"), expand = expansion(mult = c(0.02, 0))) +
  scale_y_continuous(limits = c(-100, 400), breaks = seq(-100, 400, 100)) +
  labs(x = "Survival time after resection (weeks)",
       y = "Primary tumor volume change (%)\n(Day 7 to Day 17)") +
  theme_fig

orig <- png::readPNG(file.path(WD, "data", "fig2C_ivis_image.png"))
pC <- rasterGrob(orig[205:648, 800:1420, 1:3], interpolate = TRUE)

pA <- wrap("A", "Primary tumor growth by outcome", pA)
pB <- wrap("B", "Tumor growth vs. survival (treated)", pB)
pC <- wrap("C", "IVIS one day after resection", pC)
pD <- wrap("D", "Lung signal by outcome (IVIS)", pD)
save_fig(arrangeGrob(pA, pC, pB, pD, ncol = 2, widths = c(1.1, 0.9), layout_matrix = rbind(c(1, 2), c(3, 4))),
         "Fig2.png", height = 6.2)
print(sm); print(pv); cat("Fig 2B R2 =", r2, "\n")
