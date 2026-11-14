"""Supplementary table: every tested gene set, its assigned process category,
and its survivor-vs-non-survivor enrichment at the two pre-resection timepoints.

Category assignment is a declared keyword rule, applied in the order the
categories are listed below (first match wins). The order matters: mitochondrial
rules run before translation and before stress/metabolic, so that
MITOCHONDRIAL_TRANSLATION and OXIDATIVE_PHOSPHORYLATION are classified as
mitochondrial rather than as translation or stress.

Reads (read-only):
  Manuscript_Revision_2026-07/gsea_rebuild/results/gsea_all_contrasts_long.csv
Writes:
  analysis_results/SuppTable_pathway_categories.csv
"""
import csv, collections, os, statistics

BASE = "."
GSEA = os.path.join(BASE, "Manuscript_Revision_2026-07/gsea_rebuild/results/gsea_all_contrasts_long.csv")
OUT  = os.path.join(BASE, "analysis_results/SuppTable_pathway_categories.csv")
PRE  = ["Day7_treated", "Day7_untreated", "Day18_treated", "Day18_untreated"]
ALL  = PRE + ["Week1_treated", "Week1_untreated", "Mid_treated", "Mid_untreated",
              "Late_treated", "Late_untreated"]
OUT2 = os.path.join(BASE, "analysis_results/SuppTable_figure4_catalog.csv")

RULES = [
 ("Antigen presentation & recognition", [
   "ANTIGEN","MHC","TCR_SIGNALING","T_CELL_RECEPTOR","B_CELL_RECEPTOR","BCR","CD22",
   "COSTIMULATION","PD_1","CTLA4","IMMUNOREGULATORY","CD28","IMMUNOLOGICAL_SYNAPSE"]),
 ("Inflammation & cytokine signaling", [
   "INTERFERON","INFLAMMAT","CYTOKINE","INTERLEUKIN","CHEMOKINE","TNF","NF_KB","NFKB",
   "COMPLEMENT","NEUTROPHIL","IL6","IL1","JAK_STAT","TOLL_LIKE","ANTIVIRAL","MYD88",
   "DECTIN","FCERI","FC_EPSILON","FCGR","C_TYPE_LECTIN","INFLAMMASOME","NOD_LIKE",
   "RIG_I","ACUTE_PHASE","TRAF","IRAK","IKK","RIP1","TLR","C4_AND_C2","ANTIMICROBIAL",
   "CYTOSOLIC_DNA_SENSING","STING","GPVI","PROSTANOID"]),
 ("Mitochondrial respiration", [
   "MITOCHONDRI","RESPIRATORY_ELECTRON","ELECTRON_TRANSPORT","CITRIC_ACID","TCA_CYCLE",
   "COMPLEX_I","ATP_SYNTH","OXIDATIVE_PHOSPHOR","PARKINSON","ALZHEIMER","HUNTINGTON",
   "AEROBIC","CRISTAE","PYRUVATE_METAB","BRANCHED_CHAIN"]),
 ("Translation & ribosome", [
   "TRANSLATION","RIBOSOM","RRNA","NONSENSE_MEDIATED","TRNA","EIF","SRP_DEPENDENT",
   "PEPTIDE_CHAIN","NASCENT","ELONGATION_ARREST"]),
 ("Muscle & cardiac", [
   "MUSCLE","CARDIA","CARDIO","MYOSIN","SARCOMER","DILATED","HYPERTROPHIC","STRIATED",
   "MYOGENESIS","TROPONIN"]),
 ("Stress & metabolic adaptation", [
   "HYPOXIA","AUTOPHAG","AGGREPHAGY","GLYCOLY","MTORC","UNFOLDED_PROTEIN","LYSOSOM",
   "PROTEASOME","POLYAMINE","GLYCOSYLATION","N_GLYCAN","REACTIVE_OXYGEN","OXIDATIVE_STRESS",
   "CHEMICAL_STRESS","STARVATION","HEAT_STRESS","HEAT_SHOCK","HSF1","ATTENUATION_PHASE",
   "HMOX1","CYTOPROTECT","NRF2","KEAP1","AMINO_ACID","CELLULAR_RESPONSE_TO_STRESS",
   "DETOXIF","GLUTATHIONE","BIOLOGICAL_OXIDATION","PHASE_I_FUNCTIONAL","PHASE_II_CONJUGATION",
   "CYTOCHROME_P450","UV_RESPONSE","ATF4","ENDOPLASMIC_RETICULUM_STRESS","ER_STRESS",
   "SCAVENGING_OF_HEME","IRON_UPTAKE"]),
 ("Protein homeostasis & folding", [
   "UBIQUITIN","DEUBIQUITIN","UB_SPECIFIC","SUMO","CHAPERON","PROTEIN_FOLDING","FOLDING_OF",
   "UCH_PROTEINASE","SCF_","ANAPHASE_PROMOTING","NEDDYLATION","ERAD","PREFOLDIN","TRIC_CCT",
   "CCT_TRIC","CALNEXIN","CALRETICULIN","TUBULIN_FOLDING","PROTEIN_EXPORT","QUALITY_CONTROL"]),
 ("Cancer & disease pathways", [
   "_CANCER","CARCINOMA","GLIOMA","MELANOMA","LEUKEMIA","MYELOMA","ADENOMA","NEOPLAS",
   "CYSTIC_FIBROSIS","DIABETES","LUPUS","ASTHMA","ARTHRITIS","THYROID_CANCER"]),
 ("Other immune & infection", [
   "IMMUN","LYMPH","MYELOID","PLATELET","HEMOSTASIS","COAGULATION","ALLOGRAFT","INFECTION",
   "INFLUENZA","HIV","SARS","HCMV","BACTERI","LEISHMANIA","TUBERCULOSIS","LISTERIA",
   "PHAGOCYT","GRAFT","AUTOIMMUNE","ALLERGY","VIRAL","VIRUS","MEASLES","HERPES","HEPATITIS",
   "TOXOPLASM","MALARIA","NATURAL_KILLER","HEMATOPOIETIC","NS1_MEDIATED","REV_","VPR_",
   "VIF_","NEF_","VPU_","TAT_","HOST_CELLULAR","HOST_PATHWAY"]),
 ("Cell cycle, DNA replication & death", [
   "CELL_CYCLE","MITOTIC","M_PHASE","DNA_REPAIR","DNA_REPLICATION","DNA_STRAND","SYNTHESIS_OF_DNA",
   "CHECKPOINT","G2_M","G1_S","TELOMER","CHROMATID","CHROMOSOME","SPINDLE","KINETOCHORE","E2F",
   "MYC_TARGETS","P53","APOPTO","PROGRAMMED_CELL_DEATH","NECRO","PYROPTOSIS","BH3","SENESCEN",
   "CYCLIN","DNA_DAMAGE","PRE_REPLICATIVE","ORIGIN","CDT1","CDC6","ORC_","D_LOOP",
   "HOMOLOGOUS_DNA","NUCLEAR_ENVELOPE","POLO_LIKE","AURORA"]),
 ("Transcription & RNA processing", [
   "TRANSCRIPTION","RNA_POL","SPLIC","MRNA","CAPPED_INTRON","CHROMATIN","HISTONE","EPIGEN",
   "NUCLEAR_RECEPTOR","GENE_EXPRESSION","GENE_SILENCING","RUNX","SNRNP","POLYADENYL","MIRNA",
   "ELONGATION_COMPLEX","ESTROGEN_RESPONSE","ANDROGEN"]),
 ("ECM, adhesion & migration", [
   "EXTRACELLULAR_MATRIX","ECM_","COLLAGEN","INTEGRIN","ADHESION","JUNCTION","MIGRATION",
   "EPITHELIAL_MESENCHYMAL","LAMININ","ELASTIC_FIBRE","MET_ACTIVATES","SYNDECAN","CADHERIN",
   "PROTEOGLYCAN","VASCULAR_WALL","PLATELET_ADHESION"]),
 ("Neuronal & sensory", [
   "NEURO","AXON","SYNAP","OLFACT","GABA","SEROTON","DOPAMIN","GLUTAMATE","NERVE","MYELIN",
   "OPSIN","PHOTOTRANS","TASTE","SLIT","ROBO","NETRIN","SEMA","ACETYLCHOLINE"]),
 ("Development & differentiation", [
   "DEVELOP","DIFFERENTIATION","ADIPOGENESIS","STEM_CELL","HEDGEHOG","SMO","NOTCH","WNT",
   "BETA_CATENIN","PCP_","MORPHOGEN","SPERMATO","KERATIN","CORNIFIED","EMBRYO","PLURIPOTEN",
   "ERYTHRO","MEGAKARYO","OSTEO","ANGIOGEN","EMT"]),
 ("Transport & trafficking", [
   "TRANSPORT","VESICLE","GOLGI","ENDOCYTOSIS","ENDOSOM","SECRETION","TRAFFICK","SNARE",
   "CILIUM","MEMBRANE_TRAFFIC","PROTEIN_LOCALIZATION","CLATHRIN","COPI","COPII","RAB",
   "ION_CHANNEL","ION_HOMEOSTASIS","AQUAPORIN","SLC","GLUT4","CARGO","TAIL_ANCHORED",
   "NUCLEAR_IMPORT","REABSORPTION"]),
 ("Lipid & steroid metabolism", [
   "LIPID","CHOLESTEROL","LDL","HDL","STEROID","BILE","SPHINGO","PHOSPHOLIP","TRIGLYCER",
   "PEROXISOM","EICOSANOID","PROSTAGLANDIN","FATTY_ACID","LIPOPROTEIN","KETONE","RETINOID",
   "ALDOSTERONE","RENIN_ANGIOTENSIN","PIPS","PHOSPHATIDYL","INOSITOL"]),
 ("Other metabolism", [
   "METABOLI","BIOSYNTH","CATABOL","GLUCOSE","CARBOHYDRATE","NUCLEOTIDE","PURINE","PYRIMIDINE",
   "VITAMIN","PORPHYRIN","XENOBIOTIC","DRUG_METAB","GLYCOGEN","PENTOSE","UREA","DEGRADATION",
   "GLUCONEO","SUGAR","SELENO","INSULIN","CONJUGATION","OXIDATION_BY"]),
 ("Generic signaling", [
   "SIGNALING","SIGNALLING","GPCR","RECEPTOR","KINASE","PHOSPHATASE","RHO_GTPASE","RAS_",
   "MAPK","PI3K","AKT","PTEN","RAF","CALCIUM","SECOND_MESSENGER","G_ALPHA","G_BETA","GTPASE",
   "ERBB","EGFR","FGFR","PDGF","VEGF","TGF_BETA","SMAD","HIPPO","YAP","GERANYLGERANYL",
   "TBC_RABGAP","ACTIVATION_OF"]),
]

# Categories excluded from the catalog search space, with the stated reason.
EXCLUDED = {
  "Mitochondrial respiration": "generic cellular machinery",
  "Translation & ribosome":    "generic cellular machinery",
  "Muscle & cardiac":          "host tissue at the implant site",
}


CATALOG = {}
for p in ["REACTOME_NEUTROPHIL_DEGRANULATION","HALLMARK_INFLAMMATORY_RESPONSE",
          "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
          "HALLMARK_INTERFERON_ALPHA_RESPONSE","REACTOME_INTERFERON_SIGNALING",
          "REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES",
          "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING","HALLMARK_COMPLEMENT",
          "HALLMARK_IL6_JAK_STAT3_SIGNALING","REACTOME_SIGNALING_BY_INTERLEUKINS"]:
    CATALOG[p] = "Inflammation and cytokines"
for p in ["REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION","REACTOME_TCR_SIGNALING",
          "REACTOME_DOWNSTREAM_SIGNALING_EVENTS_OF_B_CELL_RECEPTOR_BCR",
          "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
          "REACTOME_ANTIGEN_PROCESSING_UBIQUITINATION_PROTEASOME_DEGRADATION",
          "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
          "REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS",
          "REACTOME_CD22_MEDIATED_BCR_REGULATION","REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION"]:
    CATALOG[p] = "Antigen presentation and recognition"
for p in ["HALLMARK_HYPOXIA","REACTOME_AUTOPHAGY","REACTOME_METABOLISM_OF_POLYAMINES",
          "HALLMARK_MTORC1_SIGNALING",
          "HALLMARK_UNFOLDED_PROTEIN_RESPONSE","KEGG_LYSOSOME","KEGG_PROTEASOME",
          "HALLMARK_GLYCOLYSIS","REACTOME_ASPARAGINE_N_LINKED_GLYCOSYLATION",
          "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA"]:
    CATALOG[p] = "Stress and metabolic adaptation"


def category(p):
    for name, kws in RULES:
        if any(k in p for k in kws):
            return name
    return "Unclassified"


D = collections.defaultdict(dict)
meta = {}
for r in csv.DictReader(open(GSEA)):
    D[r["comparison"]][r["pathway"]] = (float(r["NES"]), float(r["padj"]))
    meta[r["pathway"]] = (r["collection"], int(r["size"]))
paths = sorted(meta)
get = lambda c, p: D[c].get(p, (float("nan"), 1.0))
nsig = lambda p: sum(1 for c in PRE if get(c, p)[1] < 0.05)
best = lambda p: min(get(c, p)[1] for c in PRE)

ranked = sorted(paths, key=best)
rank_of = {p: i + 1 for i, p in enumerate(ranked)}

os.makedirs(os.path.dirname(OUT), exist_ok=True)
nsig_all = lambda p: sum(1 for c in ALL if get(c, p)[1] < 0.05)


def contrast_cells(p):
    out = []
    for c in ALL:
        nes, q = get(c, p)
        out += ["" if nes != nes else f"{nes:.3f}", f"{q:.3g}"]
    return out


head = ["rank", "gene_set", "collection", "size", "process_category",
        "excluded_space", "exclusion_reason", "in_figure4_catalog", "catalog_category"]
for c in ALL:
    head += [f"NES_{c}", f"padj_{c}"]
head += ["n_sig_preresection_of4", "n_sig_all_of10", "best_padj_preresection"]

with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(head)
    for p in ranked:
        w.writerow([rank_of[p], p, meta[p][0], meta[p][1], category(p),
                    "Y" if category(p) in EXCLUDED else "N", EXCLUDED.get(category(p), ""),
                    "Y" if p in CATALOG else "N", CATALOG.get(p, "")]
                   + contrast_cells(p)
                   + [nsig(p), nsig_all(p), f"{best(p):.3g}"])
print("wrote", OUT, "|", len(paths), "rows")

# ---- Table Sy: the catalog only, ordered by category then rank -------------
order = {"Inflammation and cytokines": 0,
         "Antigen presentation and recognition": 1,
         "Stress and metabolic adaptation": 2}
cat_rows = sorted(CATALOG, key=lambda p: (order[CATALOG[p]], rank_of[p]))
head2 = ["catalog_category", "gene_set", "collection", "size"]
for c in ALL:
    head2 += [f"NES_{c}", f"padj_{c}"]
head2 += ["n_sig_preresection_of4", "n_sig_all_of10", "rank_of_1187"]
with open(OUT2, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(head2)
    for p in cat_rows:
        w.writerow([CATALOG[p], p, meta[p][0], meta[p][1]]
                   + contrast_cells(p)
                   + [nsig(p), nsig_all(p), rank_of[p]])
print("wrote", OUT2, "|", len(cat_rows), "rows")

# ---- summary ------------------------------------------------------------
groups = collections.defaultdict(list)
for p in paths:
    groups[category(p)].append(p)
top100 = set(ranked[:100])
print(f"\n{'category':38s} {'n':>5s} {'%univ':>6s} {'mean sig':>9s} {'%sig>=1':>8s} {'top100':>7s} {'enrich':>7s} {'catalog':>8s}")
print("-" * 100)
rows = []
for name, g in groups.items():
    rows.append((statistics.mean([nsig(p) for p in g]), name, g))
for m, name, g in sorted(rows, reverse=True):
    exp = 100 * len(g) / len(paths)
    t = sum(1 for p in g if p in top100)
    ncat = sum(1 for p in g if p in CATALOG)
    print(f"{name:38s} {len(g):5d} {exp:5.1f}% {m:9.2f} "
          f"{100*sum(1 for p in g if nsig(p)>=1)/len(g):7.0f}% {t:7d} {t/exp:6.1f}x {ncat:8d}")
print("-" * 100)
print(f"{'TOTAL':38s} {len(paths):5d}")
