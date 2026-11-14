"""Supplementary Table S4 (Supplementary Data 1): all 1,187 tested gene sets as a formatted workbook."""
import csv
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

SRC = "analysis_results/SuppTable_pathway_categories.csv"  # written by analysis_scripts/11_pathway_category_table.py
OUT = "revision_figures/Supplementary_Data_1_TableS4.xlsx"
rows = list(csv.DictReader(open(SRC)))

contrasts = [f"{tp}_{c}" for tp in ["Day7", "Day18", "Week1", "Mid", "Late"] for c in ["treated", "untreated"]]
pretty_tp = {"Day7": "Day 7", "Day18": "Day 18", "Week1": "Week 1", "Mid": "Mid", "Late": "Late"}
header = ["Rank (of 1,187)", "Gene set", "Collection", "Genes", "Process category (keyword rule)",
          "Category excluded before catalog selection", "Exclusion reason", "In top 100",
          "In 30-pathway catalog (Figs 4, 6)", "Figure category"]
for k in contrasts:
    tp, c = k.split("_")
    header += [f"NES {pretty_tp[tp]} {c}", f"padj {pretty_tp[tp]} {c}"]
header += ["Sig. (padj<0.05) in n of 4 pre-resection", "Sig. in n of 10 contrasts", "Best pre-resection padj"]

def num(x):
    try: return float(x)
    except ValueError: return x

wb = Workbook()
ws = wb.active; ws.title = "Table S4"
ws.append(header)
for r in rows:
    rank = int(r["rank"])
    line = [rank, r["gene_set"], r["collection"].title(), int(r["size"]), r["process_category"],
            "Yes" if r["excluded_space"] == "Y" else "No", r["exclusion_reason"],
            "Yes" if rank <= 100 else "No", "Yes" if r["in_figure4_catalog"] == "Y" else "No", r["catalog_category"]]
    for k in contrasts:
        line += [num(r[f"NES_{k}"]), num(r[f"padj_{k}"])]
    line += [int(r["n_sig_preresection_of4"]), int(r["n_sig_all_of10"]), num(r["best_padj_preresection"])]
    ws.append(line)

bold = Font(bold=True); fill = PatternFill("solid", fgColor="DDEBF7")
for c in ws[1]:
    c.font = bold; c.fill = fill; c.alignment = Alignment(wrap_text=True, vertical="top")
ws.row_dimensions[1].height = 45
ws.freeze_panes = "C2"
widths = {1: 9, 2: 62, 3: 11, 4: 7, 5: 30, 6: 13, 7: 26, 8: 9, 9: 13, 10: 30}
for i in range(1, len(header) + 1):
    ws.column_dimensions[get_column_letter(i)].width = widths.get(i, 11)
for row in ws.iter_rows(min_row=2, min_col=11, max_col=10 + 2 * len(contrasts)):
    for c in row:
        if isinstance(c.value, float):
            c.number_format = "0.000" if (c.column - 11) % 2 == 0 else "0.00E+00"
ws.auto_filter.ref = ws.dimensions

rd = wb.create_sheet("README", 0)
notes = [
    "Supplementary Table S4. Process category and survival enrichment for all 1,187 tested gene sets.",
    "",
    "Every Hallmark, KEGG and Reactome gene set tested by pre-ranked GSEA (fgsea on DESeq2 statistics), with its process category assigned by the keyword rule in Supplementary Table S3.",
    "Rank: position among all 1,187 gene sets ordered by the smallest adjusted p-value across the four pre-resection contrasts (Day 7 and Day 18, treated and untreated).",
    "In top 100: rank <= 100; these are the gene sets used to score process categories (Supplementary Table S3, 'In top 100').",
    "Category excluded: categories excluded before catalog selection (generic cellular machinery, etc.), with the stated reason.",
    "In 30-pathway catalog: the gene sets shown in Figures 4 and 6 (Supplementary Table S5).",
    "NES / padj: normalized enrichment score and Benjamini-Hochberg adjusted p-value for each survivor-versus-non-survivor contrast. Negative NES = enriched in survivors; positive NES = enriched in non-survivors.",
    "Mid and Late post-resection timepoints: survivors were sampled at weeks 4 and 8 after resection; non-survivors were sampled approximately 3 weeks and 1 week before death.",
]
for n in notes: rd.append([n])
rd["A1"].font = bold; rd.column_dimensions["A"].width = 160
for row in rd.iter_rows():
    for c in row: c.alignment = Alignment(wrap_text=True)
wb.save(OUT)
print("wrote", OUT, len(rows), "rows; top100 =", sum(int(r["rank"]) <= 100 for r in rows),
      "catalog =", sum(r["in_figure4_catalog"] == "Y" for r in rows))
