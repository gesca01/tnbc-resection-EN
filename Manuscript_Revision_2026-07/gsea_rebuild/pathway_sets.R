## Pathway sets and display order shared by the GSEA figure scripts.

CAT_LEVELS <- c("Inflammation and cytokines",
                "Antigen Presentation & Recognition",
                "Stress and metabolic adaptation")
CAT_COLORS <- c("Inflammation and cytokines"          = "#F4A6A6",
                "Antigen Presentation & Recognition"   = "#E15FE1",
                "Stress and metabolic adaptation"      = "#4FC978")

og_set <- data.frame(
  pathway = c(
    "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "REACTOME_INTERFERON_SIGNALING",
    "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION", "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION", "REACTOME_TCR_SIGNALING",
    "HALLMARK_HYPOXIA", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),
  module = c(rep("Inflammation and cytokines", 3),
             rep("Antigen Presentation & Recognition", 3),
             rep("Stress and metabolic adaptation", 3)),
  stringsAsFactors = FALSE)

new_set <- data.frame(
  pathway = c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "REACTOME_SIGNALING_BY_INTERLEUKINS",
    "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION", "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION", "REACTOME_TCR_SIGNALING",
    "HALLMARK_MTORC1_SIGNALING", "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),
  module = c(rep("Inflammation and cytokines", 3),
             rep("Antigen Presentation & Recognition", 3),
             rep("Stress and metabolic adaptation", 3)),
  stringsAsFactors = FALSE)

# ------------------------------------------------------------
# FULL original manuscript pathway catalog (all pathways shown in the paper's
# Figure 4A "Supporting_heatmap_curated_categories" script - not just the 9
# representative pathways used in og_set/Figure4B). Two members of new_set are
# not part of that original catalog (Signaling by Interleukins, Cellular
# Response to Hypoxia) and are appended at the end of their category so the
# "full" heatmap covers every pathway ever considered, old or new.
# highlight_pathways = the set to box/highlight (currently = new_set, i.e. the
# pathways actually in use).
#
# NOTE ON SUBSTITUTIONS: 8 pathways in the original manuscript catalog are
# KEGG_MEDICUS entries (a collection added to MSigDB after the version bundled
# with the installed msigdbr 7.5.1) or Reactome pathways since renamed --
# absent locally, which rendered as grey/blank cells. Rather than show a
# catalog that's a quarter grey, each is replaced 1:1 with a locally-available
# pathway covering the same biology (noted inline); 3 that were redundant with
# a pathway already in the catalog (Downstream-TCR ~ TCR Signaling; VCP- and
# assembly-flavored proteasome entries ~ KEGG Proteasome) are dropped rather
# than duplicated. Net: 31 pathways, all with data, zero grey cells.
# ------------------------------------------------------------
og_full_set <- data.frame(
  pathway = c(
    # Inflammation and cytokines
    "REACTOME_NEUTROPHIL_DEGRANULATION", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "REACTOME_INTERFERON_SIGNALING",
    "REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES", "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
    "HALLMARK_COMPLEMENT", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "REACTOME_SIGNALING_BY_INTERLEUKINS",
    # Antigen Presentation & Recognition
    # (Downstream-TCR-Signaling dropped - renamed/absent locally, redundant with TCR Signaling below)
    "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
    "REACTOME_TCR_SIGNALING", "REACTOME_DOWNSTREAM_SIGNALING_EVENTS_OF_B_CELL_RECEPTOR_BCR",
    "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION", "REACTOME_ANTIGEN_PROCESSING_UBIQUITINATION_PROTEASOME_DEGRADATION",
    "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
    "REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS",
    "REACTOME_CD22_MEDIATED_BCR_REGULATION", "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
    # Stress and metabolic adaptation
    "HALLMARK_HYPOXIA", "REACTOME_AUTOPHAGY",
    "REACTOME_METABOLISM_OF_POLYAMINES",          # was: Regulation of Ornithine Decarboxylase (ODC) - same polyamine-synthesis biology
    "KEGG_OXIDATIVE_PHOSPHORYLATION",              # was: KEGG_MEDICUS PINK1->Complex I electron transfer - same electron-transport-chain biology
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "KEGG_LYSOSOME",                               # was: Regulation of MITF-M genes in lysosome biogenesis/autophagy - same lysosome/autophagy biology
    "KEGG_PROTEASOME",                             # was: KEGG_MEDICUS 26S proteasome-mediated degradation - same proteasome biology
    # (VCP-to-26S-proteasome and Proteasome Assembly dropped - absent locally, redundant with KEGG Proteasome above)
    "HALLMARK_GLYCOLYSIS",                         # was: Ubiquitin-dependent degradation of Cyclin D - swapped for a core metabolic-adaptation pathway
    "REACTOME_ASPARAGINE_N_LINKED_GLYCOSYLATION",
    "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA"),
  module = c(rep("Inflammation and cytokines", 11),
             rep("Antigen Presentation & Recognition", 9),
             rep("Stress and metabolic adaptation", 11)),
  stringsAsFactors = FALSE)

highlight_pathways <- new_set$pathway
