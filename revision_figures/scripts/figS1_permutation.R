## Figure S1: label-permutation nulls for the day-18, day-7 and post-resection signatures.

source("revision_figures/scripts/common.R")
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(grid); library(gridExtra) })

SIG_ROOT  <- file.path(PROJ, "NCT_paper_pipeline/02_cache/elastic_net_day7_day18_source_match/results_EN_signatures")
POST <- file.path(PROJ, "NCT_paper_pipeline/02_cache/bmes_postresection_v2_housestyle")
vs   <- read_csv(file.path(SIG_ROOT, "validation_ci_perm_summary.csv"), show_col_types = FALSE)
post <- read_csv(file.path(POST, "summary.csv"), show_col_types = FALSE)

sig_null <- function(contrast) {
  s <- vs %>% filter(Contrast == contrast)
  list(null = read_csv(file.path(SIG_ROOT, contrast, "validation_ci_perm", "permutation_null_aucs.csv"),
                       show_col_types = FALSE)$null_auc,
       auc = s$LOOCV_AUC, p = s$Perm_p, n_perm = s$N_perms)
}
panels <- list(
  A = c(sig_null("Day18_Treated_SvNS"), title = "Resection-day signature (treated, Fig. 3B)"),
  B = c(sig_null("Day7_Treated_SvNS"),  title = "Pre-treatment signature (treated, Fig. 3C)"),
  C = list(null = read_csv(file.path(POST, "permutation_null_aucs.csv"), show_col_types = FALSE)$null_auc,
           auc = post$AUC_LOOCV, p = post$Perm_p, n_perm = post$N_perm,
           title = "Post-resection signature (Fig. 5)"))

theme_fig <- theme_classic(base_size = 11) + theme(plot.margin = margin(4, 8, 4, 8))
null_plot <- function(d) {
  ggplot(data.frame(x = d$null), aes(x)) +
    geom_histogram(bins = 25, fill = "grey78", color = "white", boundary = 0) +
    geom_vline(xintercept = d$auc, color = COL_NS, linewidth = 1.1) +
    annotate("text", x = if (median(d$null) > 0.6) 0.02 else d$auc - 0.02, y = Inf, vjust = 1.5,
             hjust = if (median(d$null) > 0.6) 0 else 1, size = 3.4, color = COL_NS,
             label = sprintf("Observed AUC %.2f\np = %.3f (%d permutations)", d$auc, d$p, d$n_perm)) +
    scale_x_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
    labs(x = "Cross-validated AUC with shuffled labels", y = "Number of permutations") + theme_fig
}
hdr <- function(letter, title) gTree(children = gList(
  textGrob(letter, x = unit(0.5, "lines"), y = 0.5, just = c("left", "center"), gp = gpar(fontsize = 16, fontface = "bold")),
  textGrob(title, x = 0.5, y = 0.5, just = "center", gp = gpar(fontsize = 12, fontface = "bold"))))
HH <- unit(1.6, "lines")
grobs <- lapply(names(panels), function(k) arrangeGrob(hdr(k, panels[[k]]$title), null_plot(panels[[k]]),
                                                      heights = unit.c(HH, unit(1, "null"))))
out <- file.path(OUT_DIR, "FigS1_permutation.png")
png(out, width = 2600, height = 950, res = 215); grid.draw(arrangeGrob(grobs = grobs, ncol = 3)); dev.off()
for (k in names(panels)) message(sprintf("%s: AUC %.3f, p %.3f, n_perm %d, null median %.2f",
                                         k, panels[[k]]$auc, panels[[k]]$p, panels[[k]]$n_perm, median(panels[[k]]$null)))
