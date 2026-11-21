## Figure S1: label-permutation nulls for the day-18, day-7 and post-resection signatures.
source("revision_figures/scripts/house_style.R")
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(grid); library(gridExtra) })

SIG_ROOT  <- file.path(PROJ, "NCT_paper_pipeline/02_cache/elastic_net_day7_day18_source_match/results_EN_signatures")
POST <- file.path(PROJ, "NCT_paper_pipeline/02_cache/bmes_postresection_v2_housestyle")
vs   <- read_csv(file.path(SIG_ROOT, "validation_ci_perm_summary.csv"), show_col_types = FALSE)
post <- read_csv(file.path(POST, "summary.csv"), show_col_types = FALSE)

v15_null <- function(contrast) {
  s <- vs %>% filter(Contrast == contrast)
  list(null = read_csv(file.path(SIG_ROOT, contrast, "validation_ci_perm", "permutation_null_aucs.csv"),
                       show_col_types = FALSE)$null_auc,
       auc = s$LOOCV_AUC, p = s$Perm_p, n_perm = s$N_perms)
}
panels <- list(
  A = c(v15_null("Day18_Treated_SvNS"), title = "Resection day (Fig. 3B)"),
  B = c(v15_null("Day7_Treated_SvNS"),  title = "Pre-treatment (Fig. 3C)"),
  C = list(null = read_csv(file.path(POST, "permutation_null_aucs.csv"), show_col_types = FALSE)$null_auc,
           auc = post$AUC_LOOCV, p = post$Perm_p, n_perm = post$N_perm,
           title = "Post-resection (Fig. 5)"))

theme_fig <- theme_house
null_plot <- function(d) {
  ggplot(data.frame(x = d$null), aes(x)) +
    geom_histogram(bins = 25, fill = "grey78", color = "white", boundary = 0) +
    geom_vline(xintercept = d$auc, color = COL_NS, linewidth = 1.1) +
    annotate("text", x = 0.02, y = Inf, vjust = 1.2, hjust = 0, size = GEOM_PT(8.5), color = COL_NS, lineheight = 0.95,
             label = sprintf("Observed AUC %.2f\np = %.3f\n(%d permutations)", d$auc, d$p, d$n_perm)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.6))) +
    scale_x_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
    labs(x = "AUC with shuffled labels", y = "Permutations") + theme_fig
}
grobs <- lapply(names(panels), function(k) wrap(k, panels[[k]]$title, null_plot(panels[[k]])))
save_fig(arrangeGrob(grobs = grobs, ncol = 3), "FigS1.png", height = 2.6)
for (k in names(panels)) message(sprintf("%s: AUC %.3f, p %.3f, n_perm %d, null median %.2f",
                                         k, panels[[k]]$auc, panels[[k]]$p, panels[[k]]$n_perm, median(panels[[k]]$null)))
