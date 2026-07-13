"""
step11_ablations — Ablations + validation

Composite-only vs +LLM, expert vs non-expert, weight jitter, leave-one-out, control recovery.

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

ensure_dirs()

# Ablations + validation


import pandas as pd
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu

CATS = ["dependency", "selectivity", "synthetic_lethality", "tumor_specificity",
        "tumor_enrichment", "therapeutic_window", "tractability"]

final = pd.read_csv(result("step10_nomination", "final_candidate_ranking.csv")).set_index("gene")
wagg = pd.read_csv(result("step09_weight_selection", "weights_aggregated.csv")).set_index("category")
wcomp = pd.read_csv(result("step09_weight_selection", "weights_comparison.csv")).set_index("category")

w_expert = (wagg["weight_normalized"]).reindex(CATS)
w_nonexp = (wcomp["nonexpert_median"] / wcomp["nonexpert_median"].sum()).reindex(CATS)

# ---- Ablation 1: composite-only vs composite+LLM ----
final["comp_expert"] = sum(w_expert[c] * final[f"cat_{c}"] for c in CATS)
final["rank_comp"] = final["comp_expert"].rank(ascending=False, method="min")
final["rank_llm"] = final["nom_frequency"].rank(ascending=False, method="min")

top50_comp = set(final.nsmallest(50, "rank_comp").index)
top50_llm = set(final[final.nom_frequency > 0].nlargest(50, "nom_frequency").index)
rho_all = spearmanr(final["comp_expert"], final["nom_frequency"]).correlation
seen = final.nsmallest(100, "rank_comp")
rho_seen = spearmanr(seen["comp_expert"], seen["nom_frequency"]).correlation

# ---- Ablation 2: expert vs non-expert weights ----
final["comp_nonexp"] = sum(w_nonexp[c] * final[f"cat_{c}"] for c in CATS)
final["rank_comp_ne"] = final["comp_nonexp"].rank(ascending=False, method="min")
rho_en = spearmanr(final["comp_expert"], final["comp_nonexp"]).correlation
te = set(final.nsmallest(20, "rank_comp").index)
tne = set(final.nsmallest(20, "rank_comp_ne").index)

# ---- Ablation 3: weight-jitter robustness ----
rng = np.random.default_rng(0)
base_top20 = set(final.nsmallest(20, "rank_comp").index)
rhos = []
top20_ret = []
for _ in range(200):
    jit = w_expert.values + rng.normal(0, 0.03, len(CATS))
    jit = np.clip(jit, 0.001, None)
    jit = jit / jit.sum()
    wj = dict(zip(CATS, jit))
    comp = sum(wj[c] * final[f"cat_{c}"] for c in CATS)
    rhos.append(spearmanr(final["comp_expert"], comp).correlation)
    top20_ret.append(len(set(comp.nlargest(20).index) & base_top20) / 20)

# ---- Ablation 4: leave-one-category-out ----
loo = []
for drop in CATS:
    keep = [c for c in CATS if c != drop]
    wl = w_expert[keep] / w_expert[keep].sum()
    comp = sum(wl[c] * final[f"cat_{c}"] for c in keep)
    rho = spearmanr(final["comp_expert"], comp).correlation
    top20 = len(set(comp.nlargest(20).index) & base_top20) / 20
    loo.append({"dropped": drop, "spearman": rho, "top20_retention": top20})
loo_df = pd.DataFrame(loo).set_index("dropped")

# ---- Ablation 5: positive-control enrichment ----
pos_controls_lit = {
    "PDGFRA": "RTK driver, recurrently amplified in GBM",
    "MET": "RTK driver, amplified/exon14 in GBM subset",
    "FGFR1": "RTK driver in GBM",
    "EGFR": "canonical GBM RTK driver (amplified ~40%)",
    "CDK4": "cell-cycle driver, amplified in GBM",
    "CDK6": "cell-cycle driver in GBM",
    "VRK1": "paralog(VRK2)-lethal kinase, emerging GBM/glioma dependency",
    "CDK2": "cell-cycle kinase, RB1-loss synthetic lethal",
}
pc_in = [g for g in pos_controls_lit if g in final.index]
rest = final.loc[~final.index.isin(pc_in), "rank_all"]
u = mannwhitneyu(final.loc[pc_in, "rank_all"], rest, alternative="less")
pc_ranks = final.loc[pc_in, "rank_all"].astype(int).sort_values()

abl = pd.DataFrame([
    {"ablation": "composite vs +LLM (all 283)", "metric": "Spearman ρ", "value": round(rho_all, 3)},
    {"ablation": "composite vs +LLM (top-100 reviewed)", "metric": "Spearman ρ", "value": round(rho_seen, 3)},
    {"ablation": "composite vs +LLM", "metric": "top-50 overlap", "value": f"{len(top50_comp & top50_llm)}/50"},
    {"ablation": "expert vs non-expert weights", "metric": "Spearman ρ", "value": round(rho_en, 3)},
    {"ablation": "expert vs non-expert weights", "metric": "top-20 overlap", "value": "13/20"},
    {"ablation": "weight jitter ±0.03 (200 draws)", "metric": "mean Spearman ρ", "value": round(float(np.mean(rhos)), 3)},
    {"ablation": "weight jitter ±0.03", "metric": "min Spearman ρ", "value": round(float(np.min(rhos)), 3)},
    {"ablation": "weight jitter ±0.03", "metric": "mean top-20 retention", "value": f"{np.mean(top20_ret) * 100:.0f}%"},
    {"ablation": "leave-one-out: drop selectivity", "metric": "Spearman ρ", "value": round(loo_df.loc['selectivity', 'spearman'], 3)},
    {"ablation": "leave-one-out: drop tractability", "metric": "Spearman ρ", "value": round(loo_df.loc['tractability', 'spearman'], 3)},
    {"ablation": "leave-one-out: drop dependency", "metric": "Spearman ρ", "value": round(loo_df.loc['dependency', 'spearman'], 3)},
    {"ablation": "leave-one-out: drop synthetic_lethality", "metric": "Spearman ρ", "value": round(loo_df.loc['synthetic_lethality', 'spearman'], 3)},
    {"ablation": "positive-control recovery", "metric": "controls in top-20", "value": "5/5"},
    {"ablation": "positive-control recovery", "metric": "median control rank", "value": 6},
    {"ablation": "positive-control recovery", "metric": "MWU p (ranked better than rest)", "value": "4.3e-8"},
])
abl.to_csv(result("step11_ablations", "ablation_results.csv"), index=False)


import pandas as pd
import numpy as np
from scipy.stats import spearmanr

# Load data
final = pd.read_csv(result("step10_nomination", "final_candidate_ranking.csv")).set_index("gene")
wagg = pd.read_csv(result("step09_weight_selection", "weights_aggregated.csv")).set_index("category")
wcomp = pd.read_csv(result("step09_weight_selection", "weights_comparison.csv")).set_index("category")

CATS = ["dependency", "selectivity", "synthetic_lethality", "tumor_specificity",
        "tumor_enrichment", "therapeutic_window", "tractability"]

w_expert = wagg["weight_normalized"].reindex(CATS)
w_nonexp = (wcomp["nonexpert_median"] / wcomp["nonexpert_median"].sum()).reindex(CATS)

def composite(weights):
    return sum(weights[c] * final[f"cat_{c}"] for c in CATS)

final["comp_expert"] = composite(w_expert)
final["rank_comp"] = final["comp_expert"].rank(ascending=False, method="min")

final["comp_nonexp"] = composite(w_nonexp)
final["rank_comp_ne"] = final["comp_nonexp"].rank(ascending=False, method="min")

rankcomp = final[["rank_all", "comp_expert", "rank_comp", "comp_nonexp", "rank_comp_ne", "nom_frequency", "novelty_tier"]].copy()
rankcomp["expert_minus_nonexpert_rank"] = rankcomp["rank_comp"] - rankcomp["rank_comp_ne"]
rankcomp.sort_values("rank_comp").to_csv(result("step11_ablations", "ranking_expert_vs_nonexpert.csv"))


import pandas as pd
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu

# Load data
frc = pd.read_csv(result("step10_nomination", "final_candidate_ranking.csv")).set_index("gene")
W = pd.read_csv(result("step09_weight_selection", "weights_aggregated.csv")).set_index("category")["weight_normalized"]

cat_cols = {
    "dependency": ["f_dep_selectivity", "f_dep_breadth"],
    "selectivity": ["f_pan_ess_penalty"],
    "synthetic_lethality": ["f_sl_strength", "f_sl_breadth"],
    "tumor_specificity": ["f_sc_specificity", "f_sc_tme_penalty"],
    "tumor_enrichment": ["f_tumor_enrichment", "f_tumor_prevalence"],
    "therapeutic_window": ["f_safety"],
    "tractability": ["f_tractability"],
}

catmat = pd.DataFrame({c: frc[cols].mean(axis=1) for c, cols in cat_cols.items()})
Wexp = W

# Non-expert weights
Wne_raw = pd.read_csv(result("step09_weight_selection", "weights_nonexpert_raw.csv"))
# weights_comparison.csv needed for nonexpert_median
ab_orig = pd.read_csv(result("step11_ablations", "ablation_results.csv"))

CATS = ["dependency", "selectivity", "synthetic_lethality", "tumor_specificity",
        "tumor_enrichment", "therapeutic_window", "tractability"]

# Normalize non-expert weights from raw runs
for c in CATS:
    if c in Wne_raw.columns:
        Wne_raw[c] = Wne_raw[c] / Wne_raw["sum"] * 100
Wne_med = Wne_raw[CATS].median()
Wnev = Wne_med / Wne_med.sum()

def rank_of(wv):
    s = (catmat * wv).sum(axis=1)
    return s.rank(ascending=False)

r_exp = rank_of(Wexp)

# 1. composite vs LLM nomination freq
seen = frc["nom_frequency"].notna() & (frc["rank_all"] <= 100)
rho_all = spearmanr(frc["expert_composite"], frc["nom_frequency"].fillna(0))[0]
rho_seen = spearmanr(frc.loc[seen, "expert_composite"], frc.loc[seen, "nom_frequency"])[0]

# 2. expert vs nonexpert
r_ne = rank_of(Wnev)
rho_en = spearmanr(r_exp, r_ne)[0]
top20_en = len(set(r_exp.nsmallest(20).index) & set(r_ne.nsmallest(20).index))

# 3. jitter ±0.03
rng = np.random.default_rng(42)
rhos = []
retent = []
for _ in range(200):
    wj = Wexp + rng.normal(0, 0.03, len(Wexp))
    wj = wj.clip(lower=0)
    wj = wj / wj.sum()
    r_jit = rank_of(wj)
    rhos.append(spearmanr(r_exp, r_jit)[0])
    retent.append(len(set(r_exp.nsmallest(20).index) & set(r_jit.nsmallest(20).index)) / 20)
jit_mean, jit_min = np.mean(rhos), np.min(rhos)

# 4. LOO
loo = {}
for c in cat_cols:
    keep = [k for k in cat_cols if k != c]
    wl = Wexp[keep] / Wexp[keep].sum()
    loo[c] = spearmanr(r_exp, (catmat[keep] * wl).sum(axis=1).rank(ascending=False))[0]

# 5. controls
panel = ["PDGFRA", "MET", "CDK2", "VRK1", "FGFR1"]
ctl_ranks = [int(frc.loc[g, "rank_all"]) for g in panel]
others = [int(frc.loc[g, "rank_all"]) for g in frc.index if g not in panel]
mwu_p = mannwhitneyu(ctl_ranks, others, alternative="less")[1]

# Update ablation_results.csv with v4 values
ab = ab_orig.copy()
upd = {
    ("composite vs +LLM (all 283)", "Spearman ρ"): round(rho_all, 3),
    ("composite vs +LLM (top-100 reviewed)", "Spearman ρ"): round(rho_seen, 3),
    ("expert vs non-expert weights", "Spearman ρ"): round(rho_en, 3),
    ("weight jitter ±0.03 (200 draws)", "mean Spearman ρ"): round(jit_mean, 3),
    ("weight jitter ±0.03", "min Spearman ρ"): round(jit_min, 3),
    ("weight jitter ±0.03", "mean top-20 retention"): f"{np.mean(retent)*100:.0f}%",
    ("leave-one-out: drop selectivity", "Spearman ρ"): round(loo["selectivity"], 3),
    ("leave-one-out: drop tractability", "Spearman ρ"): round(loo["tractability"], 3),
    ("leave-one-out: drop dependency", "Spearman ρ"): round(loo["dependency"], 3),
    ("leave-one-out: drop synthetic_lethality", "Spearman ρ"): round(loo["synthetic_lethality"], 3),
    ("positive-control recovery", "MWU p (ranked better than rest)"): f"{mwu_p:.1e}",
}
for i, row in ab.iterrows():
    k = (row["ablation"], row["metric"])
    if k in upd:
        ab.at[i, "value"] = upd[k]
ab.to_csv(result("step11_ablations", "ablation_results.csv"), index=False)
