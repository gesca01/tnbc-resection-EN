## GSVA pathway activity scores per sample (VST counts) and group-mean line plots
## across day 7 and day 18 (treated/untreated x survivor/non-survivor).

options(stringsAsFactors = FALSE)
set.seed(123)

PIPELINE_ROOT <- "./NCT_paper_pipeline"
ROOT          <- "./Manuscript_Revision_2026-07/gsea_rebuild"
.libPaths(file.path(PIPELINE_ROOT, "renv_lib"))

COUNTS_TSV   <- file.path(PIPELINE_ROOT, "01_inputs", "merged_counts.tsv")
METADATA_TSV <- file.path(PIPELINE_ROOT, "01_inputs", "merged_metadata_all.tsv")
OUT_DIR      <- file.path(ROOT, "figures", "Figure4", "pathway_activity_lineplots")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr); library(tibble)
  library(DESeq2); library(GSVA); library(msigdbr)
  library(AnnotationDbi); library(org.Mm.eg.db); library(ggplot2)
})
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

source(file.path(ROOT, "pathway_sets.R"))

strip_version <- function(x) sub("\\.\\d+$", "", x)

# ------------------------------------------------------------
# Load + map counts, same recipe as run_gsea_all_contrasts.R
# ------------------------------------------------------------
message("Loading data ...")
counts_raw <- read_tsv(COUNTS_TSV, col_types = cols())
meta_raw   <- read_tsv(METADATA_TSV, col_types = cols())
gene_col   <- names(counts_raw)[1]
counts_raw <- counts_raw %>% rename(gene_id = all_of(gene_col))

meta <- meta_raw %>%
  mutate(across(c(Timepoint, Treatment, Survival_Status, Batch), as.character)) %>%
  filter(Survival_Status %in% c("survivor", "non-survivor"),
         Timepoint %in% c("Day7", "Day18")) %>%
  mutate(Treatment = case_when(Treatment == "y" ~ "treated", Treatment == "n" ~ "untreated", TRUE ~ NA_character_),
         Survival_Status = case_when(Survival_Status == "non-survivor" ~ "non_survivor",
                                     Survival_Status == "survivor" ~ "survivor", TRUE ~ NA_character_),
         Batch = str_replace_all(Batch, "[^A-Za-z0-9_]", "_"))

sample_cols <- intersect(names(counts_raw), meta$Core_Sample_Id)
meta <- meta %>% filter(Core_Sample_Id %in% sample_cols) %>%
  mutate(Core_Sample_Id = factor(Core_Sample_Id, levels = sample_cols)) %>% arrange(Core_Sample_Id)
meta$Core_Sample_Id <- as.character(meta$Core_Sample_Id)
count_df <- counts_raw %>% select(gene_id, all_of(meta$Core_Sample_Id))

message("Mapping genes ...")
ens0 <- strip_version(count_df$gene_id)
sym <- mapIds(org.Mm.eg.db, keys = ens0, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
gk <- ifelse(!is.na(sym) & sym != "", as.character(sym), ens0)
count_df2 <- count_df %>% mutate(gene_key = gk) %>% select(gene_key, all_of(meta$Core_Sample_Id)) %>%
  group_by(gene_key) %>% summarise(across(everything(), ~ sum(as.numeric(.), na.rm = TRUE)), .groups = "drop")
cts <- as.matrix(count_df2[, -1]); rownames(cts) <- count_df2$gene_key
cts_f <- cts[rowSums(cts >= 10) >= 3, , drop = FALSE]
message("Samples: ", ncol(cts_f), "  Genes kept: ", nrow(cts_f))

# ------------------------------------------------------------
# VST-normalize (blind) across all Day7+Day18 samples together so scores
# are comparable across timepoint/treatment/survival groups.
# ------------------------------------------------------------
message("VST normalizing ...")
dds <- DESeqDataSetFromMatrix(round(cts_f), meta, design = ~ 1)
vsd <- tryCatch(vst(dds, blind = TRUE), error = function(e) varianceStabilizingTransformation(dds, blind = TRUE))
expr <- assay(vsd)

# ------------------------------------------------------------
# Gene sets restricted to the 31-pathway full catalog (og_full_set), same
# collections (H + C2:REACTOME + C2:KEGG) used to build gsea_all_contrasts_long.csv
# ------------------------------------------------------------
message("Building gene sets ...")
gs <- function(cat, sub = NULL) {
  a <- list(species = "Mus musculus", category = cat); if (!is.null(sub)) a$subcategory <- sub
  do.call(msigdbr, a) %>% select(gs_name, gene_symbol)
}
msig <- bind_rows(gs("H"), gs("C2", "CP:REACTOME"), gs("C2", "CP:KEGG"))
GDF <- msig %>% filter(gs_name %in% og_full_set$pathway) %>%
  group_by(gs_name) %>% summarise(genes = list(unique(gene_symbol)), .groups = "drop")
missing_pw <- setdiff(og_full_set$pathway, GDF$gs_name)
if (length(missing_pw) > 0) message("WARNING - pathways not found in msigdbr, skipped: ", paste(missing_pw, collapse = ", "))
GENESETS <- setNames(GDF$genes, GDF$gs_name)
message("Pathways scored: ", length(GENESETS))

# ------------------------------------------------------------
# Per-sample GSVA activity scores
# ------------------------------------------------------------
message("Running GSVA ...")
scores <- gsva(gsvaParam(expr, GENESETS, kcdf = "Gaussian", minSize = 10, maxSize = 500))
write_csv(as.data.frame(scores) %>% rownames_to_column("pathway"),
          file.path(OUT_DIR, "gsva_scores_day7_day18.csv"))

# ------------------------------------------------------------
# Group means +/- SEM: Timepoint (Pre-Tx/Resection) x Treatment x Survival
# ------------------------------------------------------------
long <- as.data.frame(scores) %>% rownames_to_column("pathway") %>%
  pivot_longer(-pathway, names_to = "Core_Sample_Id", values_to = "activity") %>%
  left_join(meta %>% select(Core_Sample_Id, Timepoint, Treatment, Survival_Status), by = "Core_Sample_Id")

group_lvls <- c("Treated / Survivor", "Treated / NonSurvivor", "Untreated / Survivor", "Untreated / NonSurvivor")
summ <- long %>%
  mutate(timepoint_label = factor(ifelse(Timepoint == "Day7", "Pre-Tx", "Resection"), levels = c("Pre-Tx", "Resection")),
         group = paste0(str_to_title(Treatment), " / ", ifelse(Survival_Status == "survivor", "Survivor", "NonSurvivor")),
         group = factor(group, levels = group_lvls)) %>%
  group_by(pathway, timepoint_label, group) %>%
  summarise(mean_activity = mean(activity), sem = sd(activity) / sqrt(n()), n = n(), .groups = "drop")

path_set <- og_full_set %>% mutate(module = factor(module, levels = CAT_LEVELS)) %>%
  filter(pathway %in% names(GENESETS))
clean_pathway <- function(x) x %>% str_replace_all("_", " ") %>% str_to_title() %>%
  str_replace("^Hallmark ", "H: ") %>% str_replace("^Reactome ", "R: ") %>% str_replace("^Kegg ", "K: ")

group_colors <- c("Treated / Survivor" = "#2166AC", "Treated / NonSurvivor" = "#B2182B",
                  "Untreated / Survivor" = "#67A9CF", "Untreated / NonSurvivor" = "#EF8A62")
group_linetype <- c("Treated / Survivor" = "solid", "Treated / NonSurvivor" = "solid",
                    "Untreated / Survivor" = "dashed", "Untreated / NonSurvivor" = "dashed")

theme_line <- theme_bw(base_size = 13) +
  theme(plot.title = element_text(size = 13.5, face = "bold", lineheight = 1.1),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 9), legend.key.width = unit(1.3, "lines"),
        panel.grid.minor = element_blank())
guide_2row <- guides(color = guide_legend(nrow = 2, byrow = TRUE),
                     linetype = guide_legend(nrow = 2, byrow = TRUE))

# ------------------------------------------------------------
# 1. One line plot per pathway
# ------------------------------------------------------------
for (i in seq_len(nrow(path_set))) {
  pw <- path_set$pathway[i]
  df <- summ %>% filter(pathway == pw)
  p <- ggplot(df, aes(timepoint_label, mean_activity, color = group, linetype = group, group = group)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.6) +
    geom_errorbar(aes(ymin = mean_activity - sem, ymax = mean_activity + sem), width = 0.08, linewidth = 0.6) +
    scale_color_manual(values = group_colors) +
    scale_linetype_manual(values = group_linetype) +
    labs(title = str_wrap(clean_pathway(pw), width = 60), x = NULL, y = "Mean GSVA activity") +
    theme_line + guide_2row
  fname <- paste0(sprintf("%02d", i), "_", str_replace_all(pw, "[^A-Za-z0-9]+", "_"), ".png")
  ggsave(file.path(OUT_DIR, fname), p, width = 6, height = 5, dpi = 300)
}
message("Wrote ", nrow(path_set), " per-pathway line plots -> ", OUT_DIR)

# ------------------------------------------------------------
# 2. Combined overview: facet_wrap by pathway, grouped/ordered by module
# ------------------------------------------------------------
summ_all <- summ %>% inner_join(path_set, by = "pathway") %>%
  mutate(pathway_label = str_wrap(clean_pathway(pathway), width = 30),
         pathway_label = factor(pathway_label, levels = unique(str_wrap(clean_pathway(path_set$pathway), width = 30))))

p_overview <- ggplot(summ_all, aes(timepoint_label, mean_activity, color = group, linetype = group, group = group)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = mean_activity - sem, ymax = mean_activity + sem), width = 0.1, linewidth = 0.4) +
  scale_color_manual(values = group_colors) +
  scale_linetype_manual(values = group_linetype) +
  facet_wrap(~ pathway_label, scales = "free_y", ncol = 5) +
  labs(title = "Pathway activity, Pre-Tx to Resection (full catalog, 31 pathways)",
       x = NULL, y = "Mean GSVA activity") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        panel.grid.minor = element_blank()) +
  guide_2row

n_row <- ceiling(nrow(path_set) / 5)
ggsave(file.path(OUT_DIR, "00_ALL_pathways_overview.png"), p_overview,
       width = 18, height = 3.2 * n_row + 1.5, dpi = 300, limitsize = FALSE)
message("Wrote combined overview -> 00_ALL_pathways_overview.png")

message("Done. Output: ", OUT_DIR)
