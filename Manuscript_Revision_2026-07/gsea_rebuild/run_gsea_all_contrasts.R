## Survivor vs non-survivor GSEA at day 7, day 18, week 1, mid and late,
## for treated and untreated mice (10 contrasts).
## DESeq2 Wald statistic -> fgsea on Hallmark, Reactome and KEGG (msigdbr 7.5.1).
## NES > 0 = higher in non-survivors.

options(stringsAsFactors = FALSE)
set.seed(123)

PIPELINE_ROOT <- "./NCT_paper_pipeline"
OUT_DIR <- "./Manuscript_Revision_2026-07/gsea_rebuild/results"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
.libPaths(file.path(PIPELINE_ROOT, "renv_lib"))

COUNTS_TSV   <- file.path(PIPELINE_ROOT, "01_inputs", "merged_counts.tsv")
METADATA_TSV <- file.path(PIPELINE_ROOT, "01_inputs", "merged_metadata_all.tsv")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr); library(tibble)
  library(DESeq2); library(fgsea); library(msigdbr); library(AnnotationDbi); library(org.Mm.eg.db)
})
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate; count <- dplyr::count

strip_version <- function(x) sub("\\.\\d+$", "", x)

reduce_design_for_subset <- function(meta_df, design_formula) {
  tt <- terms(design_formula); term_labels <- attr(tt, "term.labels")
  if (length(term_labels) == 0) return(~ 1)
  vars <- all.vars(design_formula)
  is_const <- vapply(vars, function(v) {
    if (!v %in% colnames(meta_df)) return(TRUE)
    x <- meta_df[[v]]; x <- x[!is.na(x)]; length(unique(x)) <= 1
  }, logical(1))
  const_vars <- vars[is_const]
  if (length(const_vars) == 0) return(design_formula)
  keep_terms <- term_labels[!vapply(term_labels, function(lbl)
    any(vapply(const_vars, function(cv) grepl(paste0("\\b", cv, "\\b"), lbl), logical(1))), logical(1))]
  if (length(keep_terms) == 0) return(~ 1)
  as.formula(paste("~", paste(keep_terms, collapse = " + ")))
}
build_design <- function(meta_df, design_formula) {
  d <- reduce_design_for_subset(meta_df, design_formula)
  drop_priority <- c("Batch", "Phase", "Treatment"); meta_dl <- droplevels(meta_df)
  repeat {
    tl <- attr(terms(d), "term.labels")
    if (length(tl) == 0) { d <- ~ Survival_Status; break }
    mm <- model.matrix(d, data = meta_dl)
    if (qr(mm)$rank == ncol(mm)) break
    cand <- intersect(drop_priority, tl); if (length(cand) == 0) break
    tl_new <- setdiff(tl, cand[1])
    d <- if (length(tl_new) == 0) ~ Survival_Status else as.formula(paste("~", paste(tl_new, collapse = " + ")))
    message("  design not full rank; dropped ", cand[1])
  }
  d
}
run_deseq <- function(counts, meta_df, design_formula)
  DESeq(DESeqDataSetFromMatrix(round(counts), meta_df, build_design(meta_df, design_formula)))
get_coef_name <- function(dds, pattern) { nm <- resultsNames(dds); hit <- nm[grepl(pattern, nm)]
  if (length(hit) != 1) stop("coef pattern ", pattern, " -> ", paste(nm, collapse=" | ")); hit }
mk_rank <- function(res) { r <- res$stat; names(r) <- rownames(res); sort(r[is.finite(r)], decreasing = TRUE) }
run_fgsea <- function(ranks, pathways) {
  res <- fgsea(pathways, ranks, minSize = 15, maxSize = 500)
  as_tibble(res) %>% mutate(leadingEdgePct = round(100*lengths(leadingEdge)/size,1),
                            leadingEdge = vapply(leadingEdge, function(x) paste(unlist(x), collapse=","), character(1))) %>%
    arrange(padj, desc(abs(NES))) %>% select(pathway, size, ES, NES, pval, padj, leadingEdgePct, leadingEdge)
}
subset_tc <- function(meta_full, counts_full, ...) { m <- filter(meta_full, ...); list(meta = m, counts = counts_full[, m$Core_Sample_Id, drop = FALSE]) }
make_phase <- function(meta) {
  ph <- rep(NA_character_, nrow(meta)); isS <- tolower(meta$Survival_Status)=="survivor"; isNS <- tolower(meta$Survival_Status)=="non_survivor"; tm <- meta$Timepoint
  ph[tm=="Week1"] <- "Week1"; ph[isS & tm=="Week4"] <- "Mid"; ph[isNS & tm=="Week-3"] <- "Mid"; ph[isS & tm=="Week8"] <- "Late"; ph[isNS & tm=="Week-1"] <- "Late"
  factor(ph, levels = c("Week1","Mid","Late"))
}

message("Loading data ...")
counts_raw <- read_tsv(COUNTS_TSV, col_types = cols()); meta_raw <- read_tsv(METADATA_TSV, col_types = cols())
gene_col <- names(counts_raw)[1]; counts_raw <- counts_raw %>% rename(gene_id = all_of(gene_col))
meta <- meta_raw %>% mutate(across(c(Timepoint,Treatment,Survival_Status,Batch), as.character)) %>%
  filter(Survival_Status %in% c("survivor","non-survivor")) %>%
  mutate(Treatment = case_when(Treatment=="y"~"treated", Treatment=="n"~"untreated", TRUE~NA_character_),
         Survival_Status = case_when(Survival_Status=="non-survivor"~"non_survivor", Survival_Status=="survivor"~"survivor", TRUE~NA_character_),
         Batch = str_replace_all(Batch, "[^A-Za-z0-9_]", "_"))
meta$Phase <- make_phase(meta)
meta <- meta %>% mutate(Treatment=factor(Treatment,levels=c("untreated","treated")),
                        Survival_Status=factor(Survival_Status,levels=c("survivor","non_survivor")), Batch=factor(Batch))
sample_cols <- intersect(names(counts_raw), meta$Core_Sample_Id)
meta <- meta %>% filter(Core_Sample_Id %in% sample_cols) %>% mutate(Core_Sample_Id=factor(Core_Sample_Id,levels=sample_cols)) %>% arrange(Core_Sample_Id)
meta$Core_Sample_Id <- as.character(meta$Core_Sample_Id)
count_df <- counts_raw %>% select(gene_id, all_of(meta$Core_Sample_Id))

message("Mapping genes ...")
ens0 <- strip_version(count_df$gene_id)
sym <- mapIds(org.Mm.eg.db, keys=ens0, keytype="ENSEMBL", column="SYMBOL", multiVals="first")
gk <- ifelse(!is.na(sym) & sym!="", as.character(sym), ens0)
count_df2 <- count_df %>% mutate(gene_key=gk) %>% select(gene_key, all_of(meta$Core_Sample_Id)) %>%
  group_by(gene_key) %>% summarise(across(everything(), ~sum(as.numeric(.), na.rm=TRUE)), .groups="drop")
cts <- as.matrix(count_df2[,-1]); rownames(cts) <- count_df2$gene_key
cts_f <- cts[rowSums(cts >= 10) >= 3, , drop = FALSE]
message("Genes kept: ", nrow(cts_f))

message("Building gene sets (H + Reactome + KEGG legacy) ...")
gs <- function(cat, sub=NULL) { a <- list(species="Mus musculus", category=cat); if(!is.null(sub)) a$subcategory <- sub
  do.call(msigdbr, a) %>% select(gs_name, gene_symbol) %>% mutate(collection=ifelse(is.null(sub),"HALLMARK", sub)) }
msig <- bind_rows(
  gs("H") %>% mutate(collection="HALLMARK"),
  gs("C2","CP:REACTOME") %>% mutate(collection="REACTOME"),
  gs("C2","CP:KEGG") %>% mutate(collection="KEGG"))
set_collection <- msig %>% distinct(gs_name, collection)
GDF <- msig %>% group_by(gs_name) %>% summarise(genes=list(unique(gene_symbol)), .groups="drop") %>%
  mutate(size=lengths(genes)) %>% filter(size>=15, size<=500)
GENESETS <- setNames(GDF$genes, GDF$gs_name)
message("Pathways: ", length(GENESETS))

all_long <- list(); manifest <- tibble()
run_one <- function(tag, m, c, desc) {
  if (ncol(c) < 4) { message("Skip ", tag, " (n=", ncol(c), ")"); return(invisible()) }
  m <- m %>% mutate(Treatment=droplevels(Treatment), Survival_Status=droplevels(Survival_Status), Batch=droplevels(Batch))
  if (nlevels(m$Survival_Status) < 2) { message("Skip ", tag, " (one survival class)"); return(invisible()) }
  ns <- sum(m$Survival_Status=="survivor"); nn <- sum(m$Survival_Status=="non_survivor")
  if (ns < 2 || nn < 2) { message("Skip ", tag, " (n_surv=", ns, " n_nonsurv=", nn, ")"); return(invisible()) }
  dds <- run_deseq(c, m, ~ Batch + Survival_Status)
  res <- results(dds, name = get_coef_name(dds, "^Survival_Status_"))
  g <- run_fgsea(mk_rank(res), GENESETS) %>% left_join(set_collection, by=c("pathway"="gs_name"))
  all_long[[tag]] <<- g %>% mutate(comparison=tag)
  manifest <<- bind_rows(manifest, tibble(tag=tag, n=nrow(m), n_surv=ns, n_nonsurv=nn,
                                          design=deparse(build_design(m, ~Batch+Survival_Status)),
                                          n_sig05=sum(g$padj<0.05,na.rm=TRUE), desc=desc))
  message("Done ", tag, " n=", nrow(m), " (", ns, "/", nn, ") sig05=", sum(g$padj<0.05,na.rm=TRUE))
}

for (tp in c("Day7","Day18","Week1")) for (tx in c("treated","untreated")) {
  s <- subset_tc(meta, cts_f, Timepoint==tp & Treatment==tx)
  run_one(paste0(tp,"_",tx), s$meta, s$counts, paste("Survival at", tp, tx))
}
for (ph in c("Mid","Late")) for (tx in c("treated","untreated")) {
  s <- subset_tc(meta, cts_f, !is.na(Phase) & Phase==ph & Treatment==tx)
  run_one(paste0(ph,"_",tx), s$meta, s$counts, paste("Survival at", ph, tx))
}

long <- bind_rows(all_long)
write_csv(long, file.path(OUT_DIR, "gsea_all_contrasts_long.csv"))
write_csv(manifest, file.path(OUT_DIR, "contrast_manifest.csv"))
message("\n=== MANIFEST ==="); print(as.data.frame(manifest)[,c("tag","n","n_surv","n_nonsurv","design","n_sig05")])
message("Written to ", OUT_DIR)
