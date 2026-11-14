## Fig S4: change in GSVA pathway activity from day 7 to day 18 for the 30 catalog
## pathways. Each group (treatment x outcome) is centered on its own day-7 mean.

options(stringsAsFactors = FALSE)

PROJ <- "."
PIPE <- file.path(PROJ, "NCT_paper_pipeline")
REB  <- file.path(PROJ, "Manuscript_Revision_2026-07", "gsea_rebuild")
OUT  <- file.path(PROJ, "Figures_4_6_Catalog30")
.libPaths(c(file.path(PIPE, "renv_lib"), .libPaths()))

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr); library(tibble); library(ggplot2)
})
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

source(file.path(OUT, "scripts", "pathway_sets_30.R"))

METADATA_TSV <- file.path(PIPE, "01_inputs", "merged_metadata_all.tsv")
SCORES_CSV   <- file.path(REB, "figures", "Figure4", "pathway_activity_lineplots",
                          "gsva_scores_day7_day18.csv")
if (!file.exists(SCORES_CSV)) stop("missing cached GSVA scores: ", SCORES_CSV)

## ---- metadata ----
meta <- read_tsv(METADATA_TSV, col_types = cols()) %>%
  mutate(across(c(Timepoint, Treatment, Survival_Status, Batch), as.character)) %>%
  filter(Survival_Status %in% c("survivor", "non-survivor"),
         Timepoint %in% c("Day7", "Day18")) %>%
  mutate(Treatment = case_when(Treatment == "y" ~ "treated",
                               Treatment == "n" ~ "untreated", TRUE ~ NA_character_),
         Survival_Status = case_when(Survival_Status == "non-survivor" ~ "non_survivor",
                                     Survival_Status == "survivor" ~ "survivor",
                                     TRUE ~ NA_character_)) %>%
  select(Core_Sample_Id, Timepoint, Treatment, Survival_Status)

scores <- read_csv(SCORES_CSV, show_col_types = FALSE)

group_lvls <- c("Treated / Survivor", "Treated / NonSurvivor",
                "Untreated / Survivor", "Untreated / NonSurvivor")

long <- scores %>%
  pivot_longer(-pathway, names_to = "Core_Sample_Id", values_to = "activity") %>%
  inner_join(meta, by = "Core_Sample_Id") %>%
  mutate(timepoint_label = factor(ifelse(Timepoint == "Day7", "Pre-Tx", "Resection"),
                                  levels = c("Pre-Tx", "Resection")),
         group = paste0(str_to_title(Treatment), " / ",
                        ifelse(Survival_Status == "survivor", "Survivor", "NonSurvivor")),
         group = factor(group, levels = group_lvls))

## ---- re-baseline each pathway x group to its own Pre-Tx mean ----
baseline <- long %>% filter(timepoint_label == "Pre-Tx") %>%
  group_by(pathway, group) %>% summarise(baseline_mean = mean(activity), .groups = "drop")

summ <- long %>% left_join(baseline, by = c("pathway", "group")) %>%
  mutate(activity_norm = activity - baseline_mean) %>%
  group_by(pathway, timepoint_label, group) %>%
  summarise(mean_activity = mean(activity_norm),
            sem = sd(activity_norm) / sqrt(n()), n = n(), .groups = "drop")

path_set <- catalog30 %>%
  mutate(module = factor(module, levels = CAT_LEVELS)) %>%
  filter(pathway %in% unique(scores$pathway))
stopifnot(nrow(path_set) == 30L)

clean_pathway <- function(x) x %>% str_replace_all("_", " ") %>% str_to_title() %>%
  str_replace("^Hallmark ", "H: ") %>% str_replace("^Reactome ", "R: ") %>% str_replace("^Kegg ", "K: ")

group_colors <- c("Treated / Survivor" = "#2166AC", "Treated / NonSurvivor" = "#B2182B",
                  "Untreated / Survivor" = "#67A9CF", "Untreated / NonSurvivor" = "#EF8A62")
group_linetype <- c("Treated / Survivor" = "solid", "Treated / NonSurvivor" = "solid",
                    "Untreated / Survivor" = "dashed", "Untreated / NonSurvivor" = "dashed")

summ_all <- summ %>% inner_join(path_set, by = "pathway") %>%
  mutate(pathway_label = str_wrap(clean_pathway(pathway), width = 30),
         pathway_label = factor(pathway_label,
                                levels = unique(str_wrap(clean_pathway(path_set$pathway), width = 30))))

p <- ggplot(summ_all, aes(timepoint_label, mean_activity,
                          color = group, linetype = group, group = group)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = mean_activity - sem, ymax = mean_activity + sem),
                width = 0.1, linewidth = 0.4) +
  scale_color_manual(values = group_colors) +
  scale_linetype_manual(values = group_linetype) +
  facet_wrap(~ pathway_label, scales = "free_y", ncol = 5) +
  labs(x = NULL, y = "GSVA activity change from Pre-Tx baseline") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position = "bottom", legend.title = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE))

n_row <- ceiling(nrow(path_set) / 5)
ggsave(file.path(OUT, "figures", "ActivityChange_PreTx_to_Resection_catalog30.png"),
       p, width = 18, height = 3.2 * n_row + 1.0, dpi = 300, limitsize = FALSE)
message("wrote ActivityChange_PreTx_to_Resection_catalog30.png  (",
        nrow(path_set), " pathways, ", n_row, " x 5 grid)")
