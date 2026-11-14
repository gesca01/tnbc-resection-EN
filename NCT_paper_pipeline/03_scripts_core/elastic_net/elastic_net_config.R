options(stringsAsFactors = FALSE)

PROJECT_ROOT <- "./NCT_paper_pipeline"

INPUT_COUNTS <- file.path(PROJECT_ROOT, "01_inputs", "merged_counts.tsv")
INPUT_METADATA <- file.path(PROJECT_ROOT, "01_inputs", "merged_metadata_all.tsv")

DAY7_DAY18_CACHE_DIR <- file.path(PROJECT_ROOT, "02_cache", "elastic_net_day7_day18", "results_EN_Full")
BMES_CACHE_DIR <- file.path(PROJECT_ROOT, "02_cache", "bmes_post_resection", "results_EN_Full_v2")

FIGURE_IMAGE_DIR <- file.path(PROJECT_ROOT, "04_outputs", "figure_images")

dir.create(DAY7_DAY18_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(BMES_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

