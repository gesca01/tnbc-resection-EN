## Shared paths and helpers for the figure scripts in revision_figures/.

PROJ    <- "."
WD      <- file.path(PROJ, "revision_figures")
OUT_DIR <- file.path(WD, "figures"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
TMP_DIR <- file.path(tempdir(), "wd_fig_tmp"); dir.create(TMP_DIR, showWarnings = FALSE)
.libPaths(c(file.path(PROJ, "NCT_paper_pipeline", "renv_lib"), .libPaths()))

## Evaluate a figure script up to its panel-header definition (the line defining HH),
## so its data preparation and helpers can be reused without drawing its figure.
source_upto <- function(file, env = globalenv(), keep_genes = NULL) {
  src <- readLines(file)
  src <- src[seq_len(grep("^HH <- unit\\(1.5", src)[1])]
  src <- sub("^FIG_DIR\\s*<-.*$", "FIG_DIR <- TMP_DIR", src)
  # The figure scripts re-apply a CPM>=1-in-25% filter that can drop a signature gene
  # (Col26a1 in the day-18 treated signature); always keep the signature genes.
  if (!is.null(keep_genes)) {
    assign("KEEP_GENES", keep_genes, envir = env)
    src <- sub("keep <- rowSums(cpm_mat >= 1) >= ceiling(0.25*ncol(cpm_mat))",
               "keep <- rowSums(cpm_mat >= 1) >= ceiling(0.25*ncol(cpm_mat)) | rownames(cpm_mat) %in% KEEP_GENES",
               src, fixed = TRUE)
  }
  eval(parse(text = src), envir = env)
}

COL_S <- "#1f77b4"; COL_NS <- "#d62728"
DPCR_YLAB <- "Tumor-derived tdTomato DNA\n(copies/µL)"

read_dpcr <- function() {
  d <- read.csv(file.path(WD, "data", "dpcr_tdtomato_values.csv"), stringsAsFactors = FALSE)
  d$group <- factor(d$group, levels = c("Survivor", "NonSurvivor"))
  d
}
p_stars <- function(p) if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
dpcr_test <- function(d, tp) {
  x <- d[d$timepoint == tp, ]
  suppressWarnings(wilcox.test(copies_per_uL ~ group, data = x, exact = FALSE)$p.value)
}
