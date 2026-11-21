## Figures 4, 6, S3 (pathway dot plots) and S4 (GSVA activity line plots).
source("revision_figures/scripts/house_style.R")
suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr); library(tidyr); library(tibble) })
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
CAT <- file.path(PROJ, "Figures_4_6_Catalog30")
source(file.path(CAT, "scripts", "pathway_sets_30.R"))

clean_pathway <- function(x) x %>% str_replace_all("_", " ") %>% str_to_title() %>%
  str_replace("^Hallmark ", "H: ") %>% str_replace("^Reactome ", "R: ") %>%
  str_replace("^Kegg Medicus ", "K: ") %>% str_replace("^Kegg ", "K: ")
wrap_module <- function(x) str_wrap(x, width = 14)
PATH_W <- 50

## df needs: x (column), facet_x (or NA), pathway, module, NES, padj
dotplot <- function(df, nes_title, facet_cols = TRUE) {
  lv <- str_wrap(clean_pathway(catalog30$pathway), PATH_W)
  df <- df %>% mutate(pathway_label = factor(str_wrap(clean_pathway(pathway), PATH_W), levels = rev(lv)),
                      module_label = factor(wrap_module(module), levels = wrap_module(CAT_LEVELS)),
                      sig = padj < 0.05, dot_size = pmin(-log10(pmax(padj, 1e-300)), 20))
  p <- ggplot(df, aes(x, pathway_label)) +
    geom_point(aes(size = dot_size, color = NES, fill = ifelse(sig, NES, NA_real_)),
               shape = 21, stroke = 0.7, alpha = 0.95) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "white", guide = "none") +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = nes_title,
                          guide = guide_colorbar(barheight = unit(0.6, "lines"), barwidth = unit(6, "lines"))) +
    scale_size_continuous(name = "-log10(padj)", range = c(1, 5.5), limits = c(0, 20), breaks = c(0, 5, 10, 15, 20)) +
    scale_y_discrete(position = "right") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 8, base_family = FONT) +
    theme(axis.text.x = element_text(size = PT$tick, face = "bold", color = "black"),
          axis.text.y.right = element_text(size = 8, color = "black", lineheight = 0.9),
          panel.grid.major.y = element_line(linewidth = 0.25), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", linewidth = LINE_W),
          strip.background = element_rect(fill = "grey90", color = "black", linewidth = LINE_W),
          strip.text.x = element_text(size = PT$strip, face = "bold", color = "black"),
          strip.text.y.left = element_text(size = 8, angle = 0, face = "bold", color = "black"),
          strip.placement = "outside",
          legend.title = element_text(size = 8, face = "bold", color = "black"),
          legend.text = element_text(size = 8, color = "black"),
          legend.position = "bottom", legend.box = "horizontal", panel.spacing = unit(3, "pt"))
  if (facet_cols) p + facet_grid(module_label ~ facet_x, scales = "free_y", space = "free_y", switch = "y")
  else p + facet_grid(module_label ~ ., scales = "free_y", space = "free_y", switch = "y")
}
save_gg <- function(p, name, height) save_fig(ggplotGrob(p), name, height = height)

## ---- Figures 4 and 6 ----
long <- read_csv(file.path(PROJ, "Manuscript_Revision_2026-07/gsea_rebuild/results/gsea_all_contrasts_long.csv"),
                 show_col_types = FALSE) %>% select(pathway, NES, padj, comparison)
contrast_df <- function(spec, x_levels) {
  lapply(seq_len(nrow(spec)), function(i) long %>% filter(comparison == spec$key[i]) %>%
           inner_join(catalog30, by = "pathway") %>%
           mutate(x = spec$x[i], facet_x = spec$treatment[i])) %>% bind_rows() %>%
    mutate(x = factor(x, levels = x_levels), facet_x = factor(facet_x, levels = c("Treated", "Untreated")))
}
f4 <- tribble(~x, ~treatment, ~key,
  "Pre-Tx", "Treated", "Day7_treated", "Pre-Tx", "Untreated", "Day7_untreated",
  "Resection", "Treated", "Day18_treated", "Resection", "Untreated", "Day18_untreated")
f6 <- tribble(~x, ~treatment, ~key,
  "Week 1", "Treated", "Week1_treated", "Week 1", "Untreated", "Week1_untreated",
  "Mid", "Treated", "Mid_treated", "Mid", "Untreated", "Mid_untreated",
  "Late", "Treated", "Late_treated", "Late", "Untreated", "Late_untreated")
NES_T <- "NES (blue = survivor ↑, red = non-survivor ↑)"
save_gg(dotplot(contrast_df(f4, c("Pre-Tx", "Resection")), NES_T) + theme(legend.title.position = "top"),
        "Fig4.png", height = 8.6)
save_gg(dotplot(contrast_df(f6, c("Week 1", "Mid", "Late")), NES_T) + theme(legend.title.position = "top"),
        "Fig6.png", height = 8.6)

## ---- Figure S3 (day-7 signature transfer) ----
CL <- c("EN-predicted split\n(untreated, n=11)", "Actual survival\n(untreated, n=11)", "Actual survival\n(treated, n=12)")
s3 <- read_csv(file.path(PROJ, "Supplementary_Figure_Day7_Transfer/tables/S_gsea_catalog_NES.csv"), show_col_types = FALSE) %>%
  filter(pathway %in% catalog30$pathway) %>% mutate(x = factor(contrast, levels = CL))
stopifnot(n_distinct(s3$pathway) == 30L)
save_gg(dotplot(s3, "NES (blue = (predicted) survivor \u2191,\nred = (predicted) non-survivor \u2191)", facet_cols = FALSE) +
          theme(axis.text.x = element_text(size = 8, lineheight = 0.9), legend.title.position = "top"),
        "FigS3.png", height = 8.6)

## ---- Figure S4 (GSVA activity change, pre-treatment to resection) ----
PIPE <- file.path(PROJ, "NCT_paper_pipeline")
meta <- read_tsv(file.path(PIPE, "01_inputs", "merged_metadata_all.tsv"), col_types = cols()) %>%
  mutate(across(c(Timepoint, Treatment, Survival_Status), as.character)) %>%
  filter(Survival_Status %in% c("survivor", "non-survivor"), Timepoint %in% c("Day7", "Day18")) %>%
  mutate(Treatment = ifelse(Treatment == "y", "Treated", "Untreated"),
         Outcome = ifelse(Survival_Status == "survivor", "Survivor", "Non-survivor")) %>%
  select(Core_Sample_Id, Timepoint, Treatment, Outcome)
scores <- read_csv(file.path(PROJ, "Manuscript_Revision_2026-07/gsea_rebuild/figures/Figure4/pathway_activity_lineplots/gsva_scores_day7_day18.csv"),
                   show_col_types = FALSE)
glv <- c("Treated / Survivor", "Treated / Non-survivor", "Untreated / Survivor", "Untreated / Non-survivor")
lng <- scores %>% pivot_longer(-pathway, names_to = "Core_Sample_Id", values_to = "activity") %>%
  inner_join(meta, by = "Core_Sample_Id") %>%
  mutate(tp = factor(ifelse(Timepoint == "Day7", "Pre-Tx", "Resection"), levels = c("Pre-Tx", "Resection")),
         group = factor(paste(Treatment, "/", Outcome), levels = glv))
summ <- lng %>% group_by(pathway, group) %>% mutate(act = activity - mean(activity[tp == "Pre-Tx"])) %>%
  group_by(pathway, tp, group) %>% summarise(m = mean(act), sem = sd(act) / sqrt(n()), .groups = "drop") %>%
  inner_join(catalog30, by = "pathway") %>%
  mutate(lab = factor(str_wrap(clean_pathway(pathway), 21), levels = unique(str_wrap(clean_pathway(catalog30$pathway), 21))))
stopifnot(n_distinct(summ$pathway) == 30L)
pS4 <- ggplot(summ, aes(tp, m, color = group, linetype = group, group = group)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_line(linewidth = 0.6) + geom_point(size = 1) +
  geom_errorbar(aes(ymin = m - sem, ymax = m + sem), width = 0.1, linewidth = 0.3) +
  scale_color_manual(values = setNames(c("#2166AC", "#B2182B", "#67A9CF", "#EF8A62"), glv)) +
  scale_linetype_manual(values = setNames(c("solid", "solid", "dashed", "dashed"), glv)) +
  scale_x_discrete(expand = expansion(add = 0.45)) +
  facet_wrap(~lab, scales = "free_y", ncol = 5) +
  labs(x = NULL, y = "GSVA activity change from Pre-Tx baseline") +
  theme_house + theme(axis.text = element_text(size = 6.5), strip.text = element_text(size = 6.5, lineheight = 0.9),
                      axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
                      legend.position = "bottom", legend.title = element_blank(), panel.spacing = unit(4, "pt"), strip.clip = "off") +
  guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))
save_gg(pS4, "FigS4.png", height = 9)
