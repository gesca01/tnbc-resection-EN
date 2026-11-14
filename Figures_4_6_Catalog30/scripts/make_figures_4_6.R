## Figures 4 and 6: NES heatmaps and dot plots for the 30-pathway catalog.
## NES < 0 = enriched in survivors, NES > 0 = enriched in non-survivors.

options(stringsAsFactors = FALSE)
.libPaths(c("./NCT_paper_pipeline/renv_lib", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(stringr); library(tibble); library(readr)
})

PROJ <- "."
HERE <- file.path(PROJ, "Figures_4_6_Catalog30")
GSEA <- file.path(PROJ, "Manuscript_Revision_2026-07", "gsea_rebuild",
                  "results", "gsea_all_contrasts_long.csv")
OUT  <- file.path(HERE, "figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source(file.path(HERE, "scripts", "pathway_sets_30.R"))

long <- read_csv(GSEA, show_col_types = FALSE) %>% select(pathway, NES, padj, comparison)

## ---- shared label helpers ----
clean_pathway <- function(x) x %>% str_replace_all("_", " ") %>% str_to_title() %>%
  str_replace("^Hallmark ", "H: ") %>% str_replace("^Reactome ", "R: ") %>%
  str_replace("^Kegg Medicus ", "K: ") %>% str_replace("^Kegg ", "K: ")
wrap_pathway   <- function(x, width = 56) str_wrap(clean_pathway(x), width = width)
wrap_module    <- function(x, width = 16) str_wrap(x, width = width)
truncate_label <- function(x, maxlen = 32)
  ifelse(nchar(x) > maxlen, paste0(substr(x, 1, maxlen - 3), "..."), x)
padj_to_size   <- function(p) pmin(-log10(pmax(p, 1e-300)), 20)

## row specs: which contrast goes on which row/column of each figure
f4_rows <- tribble(~timepoint_label, ~treatment,  ~key,               ~heat_label,
  "Pre-Tx",    "Treated",   "Day7_treated",    "Pre-Tx Treated",
  "Pre-Tx",    "Untreated", "Day7_untreated",  "Pre-Tx Untreated",
  "Resection", "Treated",   "Day18_treated",   "Resection Treated",
  "Resection", "Untreated", "Day18_untreated", "Resection Untreated")
f4_x <- c("Pre-Tx", "Resection")

f6_rows <- tribble(~timepoint_label, ~treatment,  ~key,               ~heat_label,
  "Week 1", "Treated",   "Week1_treated",   "Week1 Treated",
  "Week 1", "Untreated", "Week1_untreated", "Week1 Untreated",
  "Mid",    "Treated",   "Mid_treated",     "Mid Treated",
  "Mid",    "Untreated", "Mid_untreated",   "Mid Untreated",
  "Late",   "Treated",   "Late_treated",    "Late Treated",
  "Late",   "Untreated", "Late_untreated",  "Late Untreated")
f6_x <- c("Week 1", "Mid", "Late")

## ============================================================
## Panel A - full-catalog heatmap
## ============================================================
draw_heatmap <- function(rows_spec, out_hi, out_plain, width, height) {
  set_df <- catalog30 %>%
    mutate(module = factor(module, levels = CAT_LEVELS),
           pathway_label = truncate_label(clean_pathway(pathway)),
           pathway_label = factor(pathway_label, levels = unique(pathway_label)))

  df <- lapply(seq_len(nrow(rows_spec)), function(i) {
    set_df %>% mutate(row_label = rows_spec$heat_label[i]) %>%
      left_join(long %>% filter(comparison == rows_spec$key[i]) %>% select(pathway, NES, padj),
                by = "pathway")
  }) %>% bind_rows() %>%
    mutate(row_label = factor(row_label, levels = rev(rows_spec$heat_label)),
           sig  = !is.na(padj) & padj < 0.05,
           star = ifelse(sig, "*", ""),
           is_highlight = pathway %in% highlight_pathways)

  n_rows <- nrow(rows_spec)
  rng <- max(2.6, ceiling(max(abs(df$NES), na.rm = TRUE) * 10) / 10)

  base <- ggplot(df, aes(pathway_label, row_label)) +
    geom_tile(aes(fill = NES), color = "white", linewidth = 0.7) +
    geom_text(aes(label = star), size = 4.4, fontface = "bold", color = "black") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-rng, rng), na.value = "grey85", name = "NES") +
    facet_grid(. ~ module, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 11, color = "black"),
          axis.text.y = element_text(size = 13, color = "black", face = "bold"),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "grey40", linewidth = 0.6, fill = NA),
          strip.background = element_rect(fill = "grey88", color = "grey40"),
          strip.text = element_text(face = "bold", size = 12.5, color = "black"),
          panel.spacing = unit(0.3, "lines"),
          legend.position = "right",
          legend.title = element_text(size = 12, color = "black", face = "bold"),
          legend.text = element_text(size = 10.5, color = "black"),
          plot.margin = margin(12, 14, 12, 12))

  hi_df <- df %>% filter(is_highlight) %>% distinct(module, pathway_label)
  p_hi <- base + geom_tile(data = hi_df, aes(x = pathway_label, y = (n_rows + 1) / 2),
                           width = 0.94, height = n_rows, fill = NA,
                           color = "black", linewidth = 2.6, inherit.aes = FALSE)

  ggsave(out_hi,    p_hi, width = width, height = height, dpi = 300)
  ggsave(out_plain, base, width = width, height = height, dpi = 300)
  message("wrote ", basename(out_hi), " / ", basename(out_plain))
}

## ============================================================
## Panel B - dot plot, category strip on the left
## ============================================================
draw_dotplot <- function(rows_spec, x_levels, path_df, out_png, width, per_row = 0.5) {
  set_df <- path_df %>%
    mutate(module = factor(module, levels = CAT_LEVELS),
           pathway_label = wrap_pathway(pathway),
           pathway_label = factor(pathway_label, levels = unique(pathway_label)))

  df <- lapply(seq_len(nrow(rows_spec)), function(i) {
    long %>% filter(comparison == rows_spec$key[i]) %>%
      inner_join(set_df, by = "pathway") %>%
      mutate(timepoint_label = rows_spec$timepoint_label[i], treatment = rows_spec$treatment[i])
  }) %>% bind_rows() %>%
    mutate(sig = padj < 0.05, dot_size = padj_to_size(padj),
           treatment = factor(treatment, levels = c("Treated", "Untreated")),
           timepoint_label = factor(timepoint_label, levels = x_levels),
           module_label = factor(wrap_module(as.character(module)), levels = wrap_module(CAT_LEVELS)),
           pathway_label = factor(pathway_label, levels = rev(levels(set_df$pathway_label))))

  n_pathways <- nlevels(droplevels(df$pathway_label))

  p <- ggplot(df, aes(timepoint_label, pathway_label)) +
    geom_point(aes(size = dot_size, color = NES, fill = ifelse(sig, NES, NA_real_)),
               shape = 21, stroke = 1.0, alpha = 0.95) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         na.value = "white", guide = "none") +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          name = "NES\nblue = survivor ↑\nred = non-survivor ↑") +
    scale_size_continuous(name = "-log10(padj)", range = c(2.2, 10)) +
    scale_y_discrete(position = "right") +
    facet_grid(module_label ~ treatment, scales = "free_y", space = "free_y", switch = "y") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(size = 14, color = "black", face = "bold"),
          axis.text.y.right = element_text(size = 14.5, color = "black", lineheight = 0.95),
          panel.grid.major.y = element_line(linewidth = 0.3), panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey88", color = "grey40"),
          strip.text.x = element_text(size = 15, face = "bold", color = "black"),
          strip.text.y.left = element_text(size = 12.5, angle = 0, face = "bold", color = "black"),
          strip.placement = "outside",
          legend.title = element_text(size = 13, color = "black", face = "bold"),
          legend.text = element_text(size = 11.5, color = "black"),
          legend.position = "bottom", legend.box = "vertical",
          panel.spacing.y = unit(0.6, "lines"))

  ggsave(out_png, p, width = width, height = 4.2 + n_pathways * per_row, dpi = 300, limitsize = FALSE)
  message("wrote ", basename(out_png), "  (", n_pathways, " pathways)")
}

rep9 <- catalog30 %>% filter(pathway %in% highlight_pathways)

## ---- Figure 4 ----
draw_heatmap(f4_rows,
             file.path(OUT, "Figure4A_heatmap_catalog30_highlighted.png"),
             file.path(OUT, "Figure4A_heatmap_catalog30_plain.png"),
             width = 19.0, height = 5.0)
draw_dotplot(f4_rows, f4_x, catalog30,
             file.path(OUT, "Figure4B_dotplot_catalog30_groupleft.png"),
             width = 14)
draw_dotplot(f4_rows, f4_x, rep9,
             file.path(OUT, "Figure4B_dotplot_representative9.png"),
             width = 14, per_row = 0.85)

## ---- Figure 6 ----
draw_heatmap(f6_rows,
             file.path(OUT, "Figure6A_heatmap_catalog30_highlighted.png"),
             file.path(OUT, "Figure6A_heatmap_catalog30_plain.png"),
             width = 19.0, height = 6.8)
draw_dotplot(f6_rows, f6_x, catalog30,
             file.path(OUT, "Figure6B_dotplot_catalog30_groupleft.png"),
             width = 16)
draw_dotplot(f6_rows, f6_x, rep9,
             file.path(OUT, "Figure6B_dotplot_representative9.png"),
             width = 16, per_row = 0.85)

message("\nAll Figure 4 / Figure 6 panels written to ", OUT)
