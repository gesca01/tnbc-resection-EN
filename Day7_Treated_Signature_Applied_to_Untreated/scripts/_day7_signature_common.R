## Shared preprocessing for the day-7 signature transfer scripts.
## Builds the 15-gene day-7 treated signature matrix (z-scored, with the treatment-arm
## offset removed, since arm is confounded with sequencing run position at day 7)
## and the elastic-net linear predictor for each sample.

build_day7_signature_space <- function(
    proj = ".") {

  pipe <- file.path(proj, "NCT_paper_pipeline")
  sig_root <- file.path(pipe, "02_cache", "elastic_net_day7_day18_source_match",
                        "results_EN_signatures")

  strip_ver   <- function(x) sub("\\.\\d+$", "", x)
  zscore_rows <- function(m) { z <- t(scale(t(m))); z[is.na(z)] <- 0; z }

  meta_all <- readr::read_tsv(file.path(pipe, "01_inputs", "merged_metadata_all.tsv"),
                              show_col_types = FALSE)
  cts_tbl  <- readr::read_tsv(file.path(pipe, "01_inputs", "merged_counts.tsv"),
                              show_col_types = FALSE)
  gene_col <- names(cts_tbl)[1]; sample_cols <- setdiff(names(cts_tbl), gene_col)
  ens <- strip_ver(as.character(cts_tbl[[gene_col]]))
  mat_raw <- as.matrix(cts_tbl[, sample_cols]); storage.mode(mat_raw) <- "numeric"
  rownames(mat_raw) <- ens
  mat_raw <- rowsum(mat_raw, rownames(mat_raw), reorder = FALSE)

  odb <- getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
  tab <- AnnotationDbi::select(odb, keys = rownames(mat_raw), keytype = "ENSEMBL",
                               columns = "SYMBOL")
  sym <- rep(NA_character_, nrow(mat_raw))
  sym[match(tab$ENSEMBL, rownames(mat_raw))] <- tab$SYMBOL
  keep_sym <- !is.na(sym) & nzchar(sym) & !grepl("Rik$", sym) & !grepl("^(Gm|Ps)", sym)
  mat_sym <- rowsum(mat_raw[keep_sym, ], sym[keep_sym], reorder = FALSE)

  m <- meta_all %>%
    dplyr::filter(Timepoint == "Day7", Survival_Status %in% c("survivor", "non-survivor")) %>%
    dplyr::mutate(
      GroupBin = ifelse(tolower(Survival_Status) == "non-survivor", "NonSurvivor", "Survivor"),
      TxLabel  = ifelse(tolower(Treatment) == "y", "Treated", "Untreated"))
  common <- intersect(colnames(mat_sym), m$Core_Sample_Id)
  m <- m %>% dplyr::filter(Core_Sample_Id %in% common) %>%
    dplyr::arrange(match(Core_Sample_Id, common)) %>% as.data.frame()

  d <- edgeR::DGEList(mat_sym[, common, drop = FALSE])
  d <- edgeR::calcNormFactors(d, "TMM")
  keep <- rowSums(edgeR::cpm(d, log = FALSE) >= 1) >= ceiling(0.25 * ncol(d))
  d <- d[keep, , keep.lib.sizes = FALSE]
  mat <- edgeR::cpm(d, log = TRUE, prior.count = 1)
  mat <- mat[matrixStats::rowVars(mat) > 0, ]

  arm   <- factor(m$TxLabel, levels = c("Treated", "Untreated"))
  mat_z <- zscore_rows(limma::removeBatchEffect(mat, batch = arm))

  sig_tab <- readr::read_csv(file.path(sig_root, "Day7_Treated_SvNS", "selected_signature_genes.csv"),
                         show_col_types = FALSE) %>% dplyr::arrange(dplyr::desc(abs(coef)))
  sig  <- intersect(sig_tab$gene, rownames(mat_z))
  beta <- sig_tab$coef[match(sig, sig_tab$gene)]
  M    <- mat_z[sig, , drop = FALSE]

  list(M = M, meta = m, eta = as.numeric(t(M) %*% beta), beta = beta, sig = sig)
}
