# tnbc-resection-EN

Analysis code for:

> Escalona G\*, Pereles RS\*, Desai G, Kim CY, Zini C, Ross BC, Glass J, Jeruss JS, Shea LD.
> *Engineered Niche-Based Prediction of Survival and Recurrence Following Neoadjuvant Chemoimmunotherapy and Resection in Triple-Negative Breast Cancer.* (in preparation). \*co-first authors

An implantable polycaprolactone scaffold (engineered pre-metastatic niche) was sampled before treatment (day 7), at resection (day 18), and after resection in 4T1 tumor-bearing BALB/c mice, treated or not with neoadjuvant chemoimmunotherapy. This repository contains the code that turns scaffold bulk RNA-seq and digital PCR data into the elastic-net survival signatures, pathway analyses and figures in the paper.

Raw sequencing reads: NCBI GEO **[accession to be added]**. Other data (tumor measurements, survival, digital PCR): University of Michigan Deep Blue Data **[DOI to be added]**.

## Layout

The folders mirror the project as it was analyzed, so every script finds its inputs at the same relative path. **Run every script from the repository root.**

```
NCT_paper_pipeline/
  01_inputs/                       merged_counts.tsv (gene counts), merged_metadata_all.tsv (sample annotations)
  03_scripts_core/elastic_net/     elastic_net_config.R (shared paths)
  03_scripts_simplified/elastic_net/
                                   run_day7_day18_signatures.R   Fig 3 signatures (treated + untreated, day 7 / day 18)
                                   validate_ci_permutation.R                              DeLong CI + label-permutation null for the above
                                   run_postresection_signature_v2.R                           Fig 5 post-resection signature
                                   figure_postresection_signature_housestyle.R                Fig 5 panel helpers (used by revision_figures)
Figure3_signature_analysis_report/scripts/
                                   14_day18_untreated_batch1_only_rerun.R                     untreated day-18, sequencing batch 1 only
                                   16_day18_untreated_batch1_contam_filtered_15gene_rerun.R   untreated day-18 signature (Fig S2A, Table S7)
                                   17_day7_untreated_contam_filtered_15gene_rerun.R           untreated day-7 signature (Fig S2B)
                                   10_figure3_housestyle.R, 19_suppfig_S1_noPCA_day18scores.R panel helpers (used by revision_figures)
Manuscript_Revision_2026-07/gsea_rebuild/
                                   run_gsea_all_contrasts.R          DESeq2 + fgsea, 10 survivor-vs-non-survivor contrasts
                                   figure_gsea_pathway_activity_lineplots.R   GSVA single-sample pathway scores
                                   pathway_sets.R
analysis_scripts/11_pathway_category_table.py   keyword rule -> process categories (Tables S3, S4)
Figures_4_6_Catalog30/scripts/     pathway_sets_30.R (30-pathway catalog), make_figures_4_6.R (Figs 4, 6),
                                   make_activity_trajectories.R (Fig S4), make_panelF_dotplot.R (Fig S3)
Day7_Treated_Signature_Applied_to_Untreated/scripts/   01, 02 + _day7_signature_common.R   day-7 signature transfer (Fig S3 inputs)
Supplementary_Figure_Day7_Transfer/scripts/            S01, S02                              Fig S3 inputs
survival/Survival_curve_bayes.R    Bayesian Cox model, treatment vs control (Fig 1E)
revision_figures/
  data/                            dpcr_tdtomato_values.csv, fig2A_tumor_percent_change.csv,
                                   fig2B_survival_vs_tumor_change.csv, fig2C_ivis_image.png (IVIS image, panel C)
  scripts/                         fig2_rebuild.R, fig3_revised.R, fig5_revised.R, figS1_permutation.R,
                                   figS2_untreated.R, build_table_S4_xlsx.py, common.R
data_raw/dpcr/                     QuantStudio Absolute Q run summaries (resection day, week 1)
R_packages.csv                     package versions used
LICENSE                            MIT
```

## Data in this repository

**Included (inputs):**
- `NCT_paper_pipeline/01_inputs/merged_counts.tsv`: gene-level read counts for all scaffold RNA-seq samples
- `NCT_paper_pipeline/01_inputs/merged_metadata_all.tsv`: sample annotations (sample ID, sequencing batch, timepoint, treatment, date of death, survival status)
- `data_raw/dpcr/`: digital PCR run summaries (resection day, week 1)
- `revision_figures/data/`: tdTomato dPCR values, primary tumor volume changes and survival times (Fig 2), and the IVIS image for Fig 2C

**Created by the scripts (not included):**
- `NCT_paper_pipeline/02_cache/`: elastic-net signatures, cross-validated scores and permutation nulls
- `Manuscript_Revision_2026-07/gsea_rebuild/results/`: GSEA results for all contrasts
- `analysis_results/`, `Day7_Treated_Signature_Applied_to_Untreated/results/`, `Supplementary_Figure_Day7_Transfer/data/` and `tables/`: intermediate tables
- `*/figures/` and `revision_figures/Supplementary_Data_1_TableS4.xlsx`: figures and supplementary tables

Run the scripts in the order below; the figure scripts read the outputs of the earlier steps.

## Running everything (in order)

```sh
# 1. elastic-net signatures (slow: grid search with leave-one-out CV)
Rscript NCT_paper_pipeline/03_scripts_simplified/elastic_net/run_day7_day18_signatures.R
Rscript NCT_paper_pipeline/03_scripts_simplified/elastic_net/validate_ci_permutation.R
Rscript Figure3_signature_analysis_report/scripts/14_day18_untreated_batch1_only_rerun.R
Rscript Figure3_signature_analysis_report/scripts/16_day18_untreated_batch1_contam_filtered_15gene_rerun.R
Rscript Figure3_signature_analysis_report/scripts/17_day7_untreated_contam_filtered_15gene_rerun.R
Rscript NCT_paper_pipeline/03_scripts_simplified/elastic_net/run_postresection_signature_v2.R

# 2. pathway analysis
Rscript Manuscript_Revision_2026-07/gsea_rebuild/run_gsea_all_contrasts.R
Rscript Manuscript_Revision_2026-07/gsea_rebuild/figure_gsea_pathway_activity_lineplots.R
python3 analysis_scripts/11_pathway_category_table.py

# 3. day-7 signature transfer (Fig S3 inputs)
Rscript Day7_Treated_Signature_Applied_to_Untreated/scripts/01_transfer_signature.R
Rscript Day7_Treated_Signature_Applied_to_Untreated/scripts/02_gsea_predicted_groups.R
Rscript Supplementary_Figure_Day7_Transfer/scripts/S01_assemble_panel_data.R
Rscript Supplementary_Figure_Day7_Transfer/scripts/S02_gsea_alignment.R

# 4. figures and tables (fast)
Rscript Figures_4_6_Catalog30/scripts/make_figures_4_6.R
Rscript Figures_4_6_Catalog30/scripts/make_panelF_dotplot.R
Rscript Figures_4_6_Catalog30/scripts/make_activity_trajectories.R
for f in fig2_rebuild fig3_revised fig5_revised figS1_permutation figS2_untreated; do Rscript revision_figures/scripts/$f.R; done
python3 revision_figures/scripts/build_table_S4_xlsx.py
Rscript survival/Survival_curve_bayes.R
```

## Methods notes

- **Normalization:** genes with ≥1 CPM in ≥25% of samples; TMM; log2-CPM; row z-score; limma batch removal where samples span two sequencing batches.
- **Day-7 / day-18 signatures:** limma moderated-t gene ranking performed **inside** each leave-one-out fold; elastic-net grid over K genes × α (ridge-like to lasso); LOOCV AUC with DeLong CI; 200-permutation label null repeating the full pipeline.
- **Post-resection signature:** treated + untreated mice at week 1, mid and late post-resection (n = 59). Day-18-anchored frozen batch correction. Genes are ranked **once on the full cohort** (not per fold); the 500-permutation null repeats the same procedure, and its median AUC is 0.85.
- **Timepoints:** for survivors, mid/late = weeks 4/8 after resection; for non-survivors, mid/late ≈ 3 and 1 weeks before death (metadata codes `Week4`/`Week8` vs `Week-3`/`Week-1`).
- **GSEA:** DESeq2 Wald statistic ranked into fgsea on MSigDB Hallmark, KEGG (legacy) and Reactome (msigdbr 7.5.1), BH-adjusted. Negative NES = enriched in survivors.

## Requirements

R 4.5 with the Bioconductor stack; exact versions are in `R_packages.csv`. `survival/Survival_curve_bayes.R` additionally needs `brms` and `posterior` (and a Stan toolchain). Python 3 with `openpyxl` for the Table S4 workbook.
