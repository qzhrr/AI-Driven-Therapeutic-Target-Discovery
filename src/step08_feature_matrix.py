"""
step08_feature_matrix — Feature matrix

Assemble 283x12 scored features into 7 categories; normalize; category scores.

Reconstructed from the Claude Science artifact lineage of the GBM small-molecule
target-discovery project. See docs/METHODS.md for the scientific description of
this step and README.md for how paths and run modes work.

Inputs / outputs use the repository-relative helpers in common.py
(raw/ cache/ result). Run cached (default) or live (GBM_LIVE=1); see common.py.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (
    ROOT, DATA, DATA_RAW, DATA_CACHE, RESULTS, STEP_DIRS,
    LIVE, use_cache, raw, cache, result, require, ensure_dirs,
)
from adapters import host  # portable shim for host.llm / host.mcp / host.current_model
from figstyle import apply_figure_style, META_GREY

ensure_dirs()

# Feature matrix assembly + normalization + scoring


import pandas as pd
import numpy as np

# --- load each layer, authoritative source for its native features ---
dep = pd.read_csv(result("step01_dependency", "glioma_candidates.csv")).set_index("gene")
sl  = pd.read_csv(result("step03_synthetic_lethality", "sl_candidate_features.csv")).set_index("gene")
sc  = pd.read_csv(result("step04_tumor_vs_tme", "singlecell_compartment_expression.csv"), index_col=0); sc.index.name = "gene"
tw  = pd.read_csv(result("step05_therapeutic_window", "therapeutic_window_features.csv")).set_index("gene")
te  = pd.read_csv(result("step06_tumor_enrichment", "tcga_tumor_enrichment.csv"), index_col=0); te.index.name = "gene"
tr  = pd.read_csv(result("step07_tractability", "tractability_features.csv")).set_index("gene")

M = pd.DataFrame(index=dep.index); M.index.name = "gene"
# raw columns (native source)
M["glioma_mean"] = dep["glioma_mean"]; M["cohens_d"] = dep["cohens_d"]; M["dep_fdr"] = dep["fdr"]
M["glioma_dep_fraction"] = dep["glioma_dep_fraction"]; M["pan_ess_fraction"] = dep["pan_ess_fraction"]
M["is_common_essential"] = dep["is_common_essential"]; M["archetype"] = dep["archetype"]
M["sl_effect"] = sl["sl_genotype_max_effect_size"]; M["sl_n_contexts"] = sl["sl_n_contexts"]
M["sl_best_context"] = sl["sl_best_context"]; M["sl_context_type"] = sl["sl_context_type"]
M["is_SL_partner_of_lost_TSG"] = sl["is_SL_partner_of_lost_TSG"]; M["lost_TSG_partner"] = sl["lost_TSG_partner"]
M["sc_malignant_specificity"] = sc["sc_malignant_specificity"]; M["sc_tme_dominated_flag"] = sc["sc_tme_dominated_flag"]
M["sc_malignant_mean"] = sc["sc_malignant_mean"]; M["sc_measured"] = sc["sc_measured"]
M["sc_compartment"] = tw["sc_compartment"]; M["sc_specificity_tier"] = tw["sc_specificity_tier"]
M["tcga_vs_normal_log2fc"] = te["tcga_vs_normal_log2fc"]; M["tcga_pct_over_normal"] = te["tcga_pct_over_normal"]
M["tcga_gbm_median_log2tpm"] = te["tcga_gbm_median_log2tpm"]
M["safety_composite"] = tw["safety_composite"]; M["safety_vital_organ"] = tw["safety_vital_organ"]
M["safety_brain_bulk"] = tw["safety_brain_bulk"]; M["safety_brain_neuron"] = tw["safety_brain_neuron"]
M["gtex_brain_max"] = tw["gtex_brain_max"]; M["lake_neuron_max"] = tw["lake_neuron_max"]
M["cns_high_risk"] = tw["cns_high_risk"]; M["vital_high_risk"] = tw["vital_high_risk"]; M["min_subgroup_n"] = tw["min_subgroup_n"]
M["tractability_score"] = tr["tractability_score"]; M["novelty_score"] = tr["novelty_score"]
M["has_inhibitor_moa"] = tr["has_inhibitor_moa"]; M["chembl_action_types"] = tr["chembl_action_types"]
M["chembl_n_potent"] = tr["chembl_n_potent"]; M["dgidb_n_drugs"] = tr["dgidb_n_drugs"]; M["ot_drug_count"] = tr["ot_drug_count"]
M["ot_sm_approved_drug"] = tr["ot_sm_approved_drug"]
M["uniprot_kinase"] = tr["uniprot_kinase"]; M["uniprot_enzyme"] = tr["uniprot_enzyme"]
M["uniprot_membrane"] = tr["uniprot_membrane"]; M["uniprot_loc"] = tr["uniprot_loc"]

M.to_csv(result("step08_feature_matrix", "feature_matrix_raw.csv"))

def rank01(s):  # rank-percentile to [0,1], higher input -> higher output
    return s.rank(pct=True)

def mm01(s):    # min-max, zero-preserving (for SL where 0 = no SL)
    lo, hi = s.min(), s.max()
    return (s - lo) / (hi - lo) if hi > lo else s * 0

F = pd.DataFrame(index=M.index)
# ---- Dependency ----
F["f_dep_selectivity"] = rank01(-M["cohens_d"])          # more-negative d = more selective
F["f_dep_breadth"]     = M["glioma_dep_fraction"].clip(0, 1)   # already [0,1]
# ---- Selectivity (pan-essential penalty) ----
F["f_pan_ess_penalty"] = 1.0 - M["pan_ess_fraction"].clip(0, 1)  # invert: low pan-ess = good
# ---- Synthetic-lethality (bonus) ----
F["f_sl_strength"] = mm01(M["sl_effect"])                # 0 preserved for the 216 with no SL
F["f_sl_breadth"]  = mm01(M["sl_n_contexts"].astype(float))
# ---- Tumor-cell specificity ----
spec = M["sc_malignant_specificity"].copy()
spec_norm = rank01(spec)                                 # rank among measured
spec_norm[M["sc_measured"].astype(bool) == False] = 0.5    # impute neutral
F["f_sc_specificity"] = spec_norm
tme = M["sc_tme_dominated_flag"].fillna(False).astype(bool)
F["f_sc_tme_penalty"] = np.where(tme, 0.0, 1.0)            # TME-dominated = 0, else 1
F.loc[M["sc_measured"].astype(bool) == False, "f_sc_tme_penalty"] = 0.5  # unknown = neutral
# ---- Tumor enrichment ----
F["f_tumor_enrichment"] = rank01(M["tcga_vs_normal_log2fc"])
F["f_tumor_prevalence"] = M["tcga_pct_over_normal"].clip(0, 1)
# ---- Therapeutic window ----
F["f_safety"] = M["safety_composite"].clip(0, 1)
# ---- Tractability ----
F["f_tractability"] = M["tractability_score"].clip(0, 1)
# ---- Novelty ----
F["f_novelty"] = M["novelty_score"].clip(0, 1)

assert F.isna().sum().sum() == 0, F.isna().sum()[F.isna().sum() > 0]

# category map (7 categories)
CATS = {
    "dependency": ["f_dep_selectivity", "f_dep_breadth"],
    "selectivity": ["f_pan_ess_penalty"],
    "synthetic_lethality": ["f_sl_strength", "f_sl_breadth"],
    "tumor_specificity": ["f_sc_specificity", "f_sc_tme_penalty"],
    "tumor_enrichment": ["f_tumor_enrichment", "f_tumor_prevalence"],
    "therapeutic_window": ["f_safety"],
    "tractability": ["f_tractability"],
}
# category scores = mean of member features
for cat, cols in CATS.items():
    F[f"cat_{cat}"] = F[cols].mean(axis=1)
catcols = [f"cat_{c}" for c in CATS]
# equal-weight baseline composite over the 7 categories (novelty excluded from baseline)
F["baseline_composite"] = F[catcols].mean(axis=1)

F.to_csv(result("step08_feature_matrix", "feature_matrix_normalized.csv"))

# combined parquet: raw + normalized + categories + baseline
full = M.join(F, rsuffix="_norm")
full.to_parquet(result("step08_feature_matrix", "feature_matrix.parquet"))


import pandas as pd
import numpy as np

dep = pd.read_csv(result("step01_dependency", "glioma_candidates.csv")).set_index("gene")
sl  = pd.read_csv(result("step03_synthetic_lethality", "sl_candidate_features.csv")).set_index("gene")
sc  = pd.read_csv(result("step04_tumor_vs_tme", "singlecell_compartment_expression.csv"), index_col=0); sc.index.name = "gene"
tw  = pd.read_csv(result("step05_therapeutic_window", "therapeutic_window_features.csv")).set_index("gene")
te  = pd.read_csv(result("step06_tumor_enrichment", "tcga_tumor_enrichment.csv"), index_col=0); te.index.name = "gene"
tr  = pd.read_csv(result("step07_tractability", "tractability_features.csv")).set_index("gene")

M = pd.DataFrame(index=dep.index); M.index.name = "gene"
M["glioma_mean"] = dep["glioma_mean"]; M["cohens_d"] = dep["cohens_d"]; M["dep_fdr"] = dep["fdr"]
M["glioma_dep_fraction"] = dep["glioma_dep_fraction"]; M["pan_ess_fraction"] = dep["pan_ess_fraction"]
M["is_common_essential"] = dep["is_common_essential"]; M["archetype"] = dep["archetype"]
M["sl_effect"] = sl["sl_genotype_max_effect_size"]; M["sl_n_contexts"] = sl["sl_n_contexts"]
M["sl_best_context"] = sl["sl_best_context"]; M["sl_context_type"] = sl["sl_context_type"]
M["is_SL_partner_of_lost_TSG"] = sl["is_SL_partner_of_lost_TSG"]; M["lost_TSG_partner"] = sl["lost_TSG_partner"]
M["sc_malignant_specificity"] = sc["sc_malignant_specificity"]; M["sc_tme_dominated_flag"] = sc["sc_tme_dominated_flag"]
M["sc_malignant_mean"] = sc["sc_malignant_mean"]; M["sc_measured"] = sc["sc_measured"]
M["sc_compartment"] = tw["sc_compartment"]; M["sc_specificity_tier"] = tw["sc_specificity_tier"]
M["tcga_vs_normal_log2fc"] = te["tcga_vs_normal_log2fc"]; M["tcga_pct_over_normal"] = te["tcga_pct_over_normal"]
M["tcga_gbm_median_log2tpm"] = te["tcga_gbm_median_log2tpm"]
M["safety_composite"] = tw["safety_composite"]; M["safety_vital_organ"] = tw["safety_vital_organ"]
M["safety_brain_bulk"] = tw["safety_brain_bulk"]; M["safety_brain_neuron"] = tw["safety_brain_neuron"]
M["gtex_brain_max"] = tw["gtex_brain_max"]; M["lake_neuron_max"] = tw["lake_neuron_max"]
M["cns_high_risk"] = tw["cns_high_risk"]; M["vital_high_risk"] = tw["vital_high_risk"]; M["min_subgroup_n"] = tw["min_subgroup_n"]
M["tractability_score"] = tr["tractability_score"]; M["novelty_score"] = tr["novelty_score"]
M["has_inhibitor_moa"] = tr["has_inhibitor_moa"]; M["chembl_action_types"] = tr["chembl_action_types"]
M["chembl_n_potent"] = tr["chembl_n_potent"]; M["dgidb_n_drugs"] = tr["dgidb_n_drugs"]; M["ot_drug_count"] = tr["ot_drug_count"]
M["ot_sm_approved_drug"] = tr["ot_sm_approved_drug"]
M["uniprot_kinase"] = tr["uniprot_kinase"]; M["uniprot_enzyme"] = tr["uniprot_enzyme"]
M["uniprot_membrane"] = tr["uniprot_membrane"]; M["uniprot_loc"] = tr["uniprot_loc"]

F = pd.DataFrame(index=M.index)
F["f_dep_selectivity"] = rank01(-M["cohens_d"])
F["f_dep_breadth"] = M["glioma_dep_fraction"].clip(0, 1)
F["f_pan_ess_penalty"] = 1.0 - M["pan_ess_fraction"].clip(0, 1)
F["f_sl_strength"] = mm01(M["sl_effect"])
F["f_sl_breadth"] = mm01(M["sl_n_contexts"].astype(float))
spec = M["sc_malignant_specificity"].copy()
spec_norm = rank01(spec)
spec_norm[M["sc_measured"].astype(bool) == False] = 0.5
F["f_sc_specificity"] = spec_norm
tme = M["sc_tme_dominated_flag"].fillna(False).astype(bool)
F["f_sc_tme_penalty"] = np.where(tme, 0.0, 1.0)
F.loc[M["sc_measured"].astype(bool) == False, "f_sc_tme_penalty"] = 0.5
F["f_tumor_enrichment"] = rank01(M["tcga_vs_normal_log2fc"])
F["f_tumor_prevalence"] = M["tcga_pct_over_normal"].clip(0, 1)
F["f_safety"] = M["safety_composite"].clip(0, 1)
F["f_tractability"] = M["tractability_score"].clip(0, 1)
F["f_novelty"] = M["novelty_score"].clip(0, 1)

CATS = {
    "dependency": ["f_dep_selectivity", "f_dep_breadth"],
    "selectivity": ["f_pan_ess_penalty"],
    "synthetic_lethality": ["f_sl_strength", "f_sl_breadth"],
    "tumor_specificity": ["f_sc_specificity", "f_sc_tme_penalty"],
    "tumor_enrichment": ["f_tumor_enrichment", "f_tumor_prevalence"],
    "therapeutic_window": ["f_safety"],
    "tractability": ["f_tractability"],
}

for cat, cols in CATS.items():
    F[f"cat_{cat}"] = F[cols].mean(axis=1)

catcols = [f"cat_{c}" for c in CATS]
F["baseline_composite"] = F[catcols].mean(axis=1)

scored = ["f_dep_selectivity", "f_dep_breadth", "f_pan_ess_penalty", "f_sl_strength", "f_sl_breadth",
          "f_sc_specificity", "f_sc_tme_penalty", "f_tumor_enrichment", "f_tumor_prevalence",
          "f_safety", "f_tractability", "f_novelty"]
cats = [c for c in F.columns if c.startswith("cat_")]
F[scored + cats + ["baseline_composite"]].round(4).to_csv(result("step08_feature_matrix", "feature_matrix_scored.csv"))


import pandas as pd
import numpy as np
import json

# Global state
apply_figure_style()

# Load each layer
dep = pd.read_csv(result("step01_dependency", "glioma_candidates.csv")).set_index("gene")
sl  = pd.read_csv(result("step03_synthetic_lethality", "sl_candidate_features.csv")).set_index("gene")
sc  = pd.read_csv(result("step04_tumor_vs_tme", "singlecell_compartment_expression.csv"), index_col=0); sc.index.name = "gene"
tw  = pd.read_csv(result("step05_therapeutic_window", "therapeutic_window_features.csv")).set_index("gene")
te  = pd.read_csv(result("step06_tumor_enrichment", "tcga_tumor_enrichment.csv"), index_col=0); te.index.name = "gene"
tr  = pd.read_csv(result("step07_tractability", "tractability_features.csv")).set_index("gene")

M = pd.DataFrame(index=dep.index); M.index.name = "gene"
M["glioma_mean"] = dep["glioma_mean"]; M["cohens_d"] = dep["cohens_d"]; M["dep_fdr"] = dep["fdr"]
M["glioma_dep_fraction"] = dep["glioma_dep_fraction"]; M["pan_ess_fraction"] = dep["pan_ess_fraction"]
M["is_common_essential"] = dep["is_common_essential"]; M["archetype"] = dep["archetype"]
M["sl_effect"] = sl["sl_genotype_max_effect_size"]; M["sl_n_contexts"] = sl["sl_n_contexts"]
M["sl_best_context"] = sl["sl_best_context"]; M["sl_context_type"] = sl["sl_context_type"]
M["is_SL_partner_of_lost_TSG"] = sl["is_SL_partner_of_lost_TSG"]; M["lost_TSG_partner"] = sl["lost_TSG_partner"]
M["sc_malignant_specificity"] = sc["sc_malignant_specificity"]; M["sc_tme_dominated_flag"] = sc["sc_tme_dominated_flag"]
M["sc_malignant_mean"] = sc["sc_malignant_mean"]; M["sc_measured"] = sc["sc_measured"]
M["sc_compartment"] = tw["sc_compartment"]; M["sc_specificity_tier"] = tw["sc_specificity_tier"]
M["tcga_vs_normal_log2fc"] = te["tcga_vs_normal_log2fc"]; M["tcga_pct_over_normal"] = te["tcga_pct_over_normal"]
M["tcga_gbm_median_log2tpm"] = te["tcga_gbm_median_log2tpm"]
M["safety_composite"] = tw["safety_composite"]; M["safety_vital_organ"] = tw["safety_vital_organ"]
M["safety_brain_bulk"] = tw["safety_brain_bulk"]; M["safety_brain_neuron"] = tw["safety_brain_neuron"]
M["gtex_brain_max"] = tw["gtex_brain_max"]; M["lake_neuron_max"] = tw["lake_neuron_max"]
M["cns_high_risk"] = tw["cns_high_risk"]; M["vital_high_risk"] = tw["vital_high_risk"]; M["min_subgroup_n"] = tw["min_subgroup_n"]
M["tractability_score"] = tr["tractability_score"]; M["novelty_score"] = tr["novelty_score"]
M["has_inhibitor_moa"] = tr["has_inhibitor_moa"]; M["chembl_action_types"] = tr["chembl_action_types"]
M["chembl_n_potent"] = tr["chembl_n_potent"]; M["dgidb_n_drugs"] = tr["dgidb_n_drugs"]; M["ot_drug_count"] = tr["ot_drug_count"]
M["ot_sm_approved_drug"] = tr["ot_sm_approved_drug"]
M["uniprot_kinase"] = tr["uniprot_kinase"]; M["uniprot_enzyme"] = tr["uniprot_enzyme"]
M["uniprot_membrane"] = tr["uniprot_membrane"]; M["uniprot_loc"] = tr["uniprot_loc"]

M.to_csv(result("step08_feature_matrix", "feature_matrix_raw.csv"))


F = pd.DataFrame(index=M.index)
# ---- Dependency ----
F["f_dep_selectivity"] = rank01(-M["cohens_d"])          # more-negative d = more selective
F["f_dep_breadth"]     = M["glioma_dep_fraction"].clip(0, 1)   # already [0,1]
# ---- Selectivity (pan-essential penalty) ----
F["f_pan_ess_penalty"] = 1.0 - M["pan_ess_fraction"].clip(0, 1)  # invert: low pan-ess = good
# ---- Synthetic-lethality (bonus) ----
F["f_sl_strength"] = mm01(M["sl_effect"])                # 0 preserved for the 216 with no SL
F["f_sl_breadth"]  = mm01(M["sl_n_contexts"].astype(float))
# ---- Tumor-cell specificity ----
spec = M["sc_malignant_specificity"].copy()
spec_norm = rank01(spec)                                 # rank among measured
spec_norm[M["sc_measured"].astype(bool) == False] = 0.5    # impute neutral
F["f_sc_specificity"] = spec_norm
tme = M["sc_tme_dominated_flag"].fillna(False).astype(bool)
F["f_sc_tme_penalty"] = np.where(tme, 0.0, 1.0)            # TME-dominated = 0, else 1
F.loc[M["sc_measured"].astype(bool) == False, "f_sc_tme_penalty"] = 0.5  # unknown = neutral
# ---- Tumor enrichment ----
F["f_tumor_enrichment"] = rank01(M["tcga_vs_normal_log2fc"])
F["f_tumor_prevalence"] = M["tcga_pct_over_normal"].clip(0, 1)
# ---- Therapeutic window ----
F["f_safety"] = M["safety_composite"].clip(0, 1)
# ---- Tractability ----
F["f_tractability"] = M["tractability_score"].clip(0, 1)
# ---- Novelty ----
F["f_novelty"] = M["novelty_score"].clip(0, 1)

assert F.isna().sum().sum() == 0, F.isna().sum()[F.isna().sum() > 0]

# category map (7 categories)
CATS = {
    "dependency": ["f_dep_selectivity", "f_dep_breadth"],
    "selectivity": ["f_pan_ess_penalty"],
    "synthetic_lethality": ["f_sl_strength", "f_sl_breadth"],
    "tumor_specificity": ["f_sc_specificity", "f_sc_tme_penalty"],
    "tumor_enrichment": ["f_tumor_enrichment", "f_tumor_prevalence"],
    "therapeutic_window": ["f_safety"],
    "tractability": ["f_tractability"],
}
# category scores = mean of member features
for cat, cols in CATS.items():
    F[f"cat_{cat}"] = F[cols].mean(axis=1)
catcols = [f"cat_{c}" for c in CATS]
# equal-weight baseline composite over the 7 categories (novelty excluded from baseline)
F["baseline_composite"] = F[catcols].mean(axis=1)
json.dump({k: v for k, v in CATS.items()}, open(cache("categories.json"), "w"))

F.to_csv(result("step08_feature_matrix", "feature_matrix_normalized.csv"))


import pandas as pd
import numpy as np

# Load inputs
fmx = pd.read_parquet(result("step08_feature_matrix", "feature_matrix.parquet"))
if 'gene' in fmx.columns:
    fmx = fmx.set_index('gene')
fmx.index.name = 'gene'

w_agg = pd.read_csv(result("step09_weight_selection", "weights_aggregated.csv"))

cats = ['dependency', 'selectivity', 'synthetic_lethality', 'tumor_specificity',
        'tumor_enrichment', 'therapeutic_window', 'tractability']
X = fmx[[f"cat_{c}" for c in cats]].copy()
X.columns = cats
wmap = dict(zip(w_agg['category'], w_agg['weight_normalized']))
w = np.array([wmap[c] for c in cats])

# Expert composite
fmx['expert_composite'] = X.values @ w

# Fix novelty labeling
appr = fmx['ot_sm_approved_drug'] == True
sm_chem = (fmx['has_inhibitor_moa'] == True) | (fmx['chembl_n_potent'] >= 10)
tier = pd.Series('novel', index=fmx.index)
tier[sm_chem & ~appr] = 'has_sm_chemistry'
tier[appr] = 'approved_sm_drug'
fmx['novelty_tier'] = tier
fmx['novelty_corrected'] = tier.map({'novel': 1.0, 'has_sm_chemistry': 0.5, 'approved_sm_drug': 0.0})

# Positive-control anchors
POS = ['VRK1', 'MET', 'PTPN11', 'PDGFRA', 'FGFR1']
fmx['is_dependency_poscontrol'] = fmx.index.isin(POS)

# Rank
fmx = fmx.sort_values('expert_composite', ascending=False)
fmx['rank_full'] = np.arange(1, len(fmx) + 1)
gated = fmx[fmx['novelty_tier'] != 'approved_sm_drug'].copy()
gated['rank_gated'] = np.arange(1, len(gated) + 1)
gated['in_top100_shortlist'] = gated['rank_gated'] <= 100

# Save corrected master matrix
fmx.reset_index().to_parquet(result("step08_feature_matrix", "feature_matrix.parquet"), index=False)
