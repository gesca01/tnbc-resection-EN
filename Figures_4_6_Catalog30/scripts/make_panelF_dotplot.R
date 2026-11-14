## Fig S3 dot plot: GSEA of the day-7 signature transfer on the 30-pathway catalog.
## Red = up in (predicted) non-survivors; filled points = padj < 0.05.

options(stringsAsFactors = FALSE)
.libPaths(c("./NCT_paper_pipeline/renv_lib", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(stringr)
})

PROJ <- "."
SUPP <- file.path(PROJ, "Supplementary_Figure_Day7_Transfer")
OUT  <- file.path(PROJ, "Figures_4_6_Catalog30")
source(file.path(OUT, "scripts", "pathway_sets_30.R"))

gsea <- read_csv(file.path(SUPP, "tables", "S_gsea_catalog_NES.csv"), show_col_types = FALSE)

clean_pathway <- function(x) x %>% str_replace_all("_", " ") %>% str_to_title() %>%
  str_replace("^Hallmark ", "H: ") %>% str_replace("^Reactome ", "R: ") %>% str_replace("^Kegg ", "K: ")
wrap_module <- function(x) str_wrap(x, width = 16)   # matches make_figures_4_6.R

CONTRAST_LEVELS <- c("EN-predicted split\n(untreated, n=11)",
                     "Actual survival\n(untreated, n=11)",
                     "Actual survival\n(treated, n=12)")

path_levels <- str_wrap(clean_pathway(catalog30$pathway), 56)  # matches make_figures_4_6.R

g <- gsea %>%
  filter(pathway %in% catalog30$pathway) %>%
  mutate(pathway_label = factor(str_wrap(clean_pathway(pathway), 56), levels = rev(path_levels)),
         module_label  = factor(wrap_module(module), levels = wrap_module(CAT_LEVELS)),
         contrast      = factor(contrast, levels = CONTRAST_LEVELS),
         sig           = padj < 0.05,
         dot_size      = pmin(-log10(pmax(padj, 1e-300)), 20))

stopifnot(n_distinct(g$pathway) == 30L)

## Layout matches draw_dotplot() in make_figures_4_6.R: switch = "y" puts the
## category strip on the left, scale_y_discrete(position = "right") moves the
## pathway names to the right axis, strip.placement = "outside" keeps the strip
## outboard of the panel.
p <- ggplot(g, aes(contrast, pathway_label)) +
  geom_point(aes(size = dot_size, colour = NES, fill = ifelse(sig, NES, NA_real_)),
             shape = 21, stroke = 1.0, alpha = 0.95) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       na.value = "white", guide = "none") +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         name = "NES\nred = up in (predicted) non-survivor ↑\nblue = up in (predicted) survivor ↑") +
  scale_size_continuous(name = "-log10(padj)", range = c(2.2, 10)) +
  scale_y_discrete(position = "right") +
  facet_grid(module_label ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(size = 14, color = "black", face = "bold", lineheight = 0.95),
        axis.text.y.right = element_text(size = 14.5, color = "black", lineheight = 0.95),
        panel.grid.major.y = element_line(linewidth = 0.3), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey88", color = "grey40"),
        strip.text.y.left = element_text(size = 12.5, angle = 0, face = "bold", color = "black"),
        strip.placement = "outside",
        legend.title = element_text(size = 13, color = "black", face = "bold"),
        legend.text = element_text(size = 11.5, color = "black"),
        legend.position = "bottom", legend.box = "vertical",
        panel.spacing.y = unit(0.6, "lines"))

n_pathways <- nlevels(droplevels(g$pathway_label))
ggsave(file.path(OUT, "figures", "PanelF_Day7_transfer_dotplot_catalog30.png"),
       p, width = 14, height = 4.2 + n_pathways * 0.5, dpi = 300, limitsize = FALSE)
## PNG only, like the other panels in this folder. No PDF: neither the default
## pdf device nor cairo_pdf (X11 unavailable on this machine) can render the
## U+2191 arrow in the legend title, and they fail by dropping the glyph
## silently rather than erroring.
message("wrote PanelF_Day7_transfer_dotplot_catalog30.png  (30 pathways)")
