"""
step10_nomination — LLM target nomination + final candidate ranking

Two parts:

  1. Nomination. An LLM (GBM small-molecule target-discovery framing) reviews the
     top-100 candidates by expert-weighted composite, each with a compact evidence
     dossier, and nominates the ~15 strongest 100 times. Nomination frequency across
     runs is an ensemble judgment layered on top of the quantitative composite.

  2. Final ranking. Assembles the master 283-gene table: the expert-weighted
     composite (with synthetic-lethality augmented by the step03 expression-paralog
     scan), nomination frequency, and the step11a gene-masked ablation columns.

Iterative provenance (see docs/METHODS.md). In the original interactive project the
ranking was refined in three passes: an initial ranking (v1), then re-scored after
the gene-masked ablation (step11a, v3), then after folding in the expression-context
SL scan (step03, v4). The committed feature_matrix.parquet carries both the v1
provenance columns and the final augmented values, so this script recomputes the
final (v4) composite deterministically and reproduces the published ranking exactly.

Outputs (results/step10_nomination/):
  nomination_base_p2soft.csv, nomination_runs_raw.csv, nomination_cautions.csv,
  nomination_frequency.csv, final_candidate_ranking.csv

Run modes (see common.py):
  * cached (default): committed nomination runs are used as-is — no LLM calls.
  * live (GBM_LIVE=1): re-query the LLM for nomination (needs 'anthropic' +
    ANTHROPIC_API_KEY); stochastic, so it will not reproduce the committed runs exactly.
"""
import os
import sys
import json

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import LIVE, use_cache, result, ensure_dirs  # noqa: E402
from adapters import host  # noqa: E402

ensure_dirs()

STEP = "step10_nomination"
CATS = ["dependency", "selectivity", "synthetic_lethality", "tumor_specificity",
        "tumor_enrichment", "therapeutic_window", "tractability"]
CAT_COLS = {
    "dependency": ["f_dep_selectivity", "f_dep_breadth"],
    "selectivity": ["f_pan_ess_penalty"],
    "synthetic_lethality": ["f_sl_strength", "f_sl_breadth"],
    "tumor_specificity": ["f_sc_specificity", "f_sc_tme_penalty"],
    "tumor_enrichment": ["f_tumor_enrichment", "f_tumor_prevalence"],
    "therapeutic_window": ["f_safety"],
    "tractability": ["f_tractability"],
}

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
fm = pd.read_parquet(result("step08_feature_matrix", "feature_matrix.parquet"))
if "gene" in fm.columns:
    fm = fm.set_index("gene")
tr = pd.read_csv(result("step07_tractability", "tractability_features.csv")).set_index("gene")
W = pd.read_csv(result("step09_weight_selection", "weights_aggregated.csv")).set_index("category")["weight_normalized"]

# ---------------------------------------------------------------------------
# Synthetic-lethality augmentation (step03 expression-paralog scan) -> v4 composite
# ---------------------------------------------------------------------------
# Recompute the augmented SL features from the v1 base columns + step03's expression
# scan, so the composite is reproducible from committed inputs (not just inherited).
scan = pd.read_csv(result("step03_synthetic_lethality", "expression_SL_scan_candidates.csv")).set_index("gene")


def _is_sig(v):
    return (v is True) or (str(v) in ("True", "1", "1.0"))


def _mm01(s):
    s = s.astype(float)
    lo, hi = s.min(), s.max()
    return (s - lo) / (hi - lo) if hi > lo else s * 0


# All 10 name-root pairs are Ensembl-confirmed true paralogs; keep the significant ones.
ENSEMBL_PARALOGS = ["FAM50A", "CDS2", "RPP25L", "EAF1", "PTK2", "ELMO2", "CHMP3", "GOLT1B", "VRK1", "ACTL6A"]

if {"sl_effect_v1", "sl_n_contexts_v1"}.issubset(fm.columns):
    cred = {g: abs(float(scan.loc[g, "cohens_d"]))
            for g in ENSEMBL_PARALOGS if g in scan.index and _is_sig(scan.loc[g, "sig"])}
    cap = float(fm["sl_effect_v1"].max())  # winsorize at the strongest genotype SL
    sl_effect = pd.Series([max(fm.loc[g, "sl_effect_v1"], cred.get(g, 0.0)) for g in fm.index], index=fm.index)
    added = pd.Series({g: int(cred.get(g, 0.0) > fm.loc[g, "sl_effect_v1"]) for g in fm.index})
    fm["f_sl_strength"] = _mm01(sl_effect.clip(upper=cap))
    fm["f_sl_breadth"] = _mm01((fm["sl_n_contexts_v1"] + added).astype(float))
    fm["sl_effect"] = sl_effect
    fm["sl_n_contexts"] = fm["sl_n_contexts_v1"] + added
    fm["sl_expr_partner"] = [scan.loc[g, "top_expr_partner"] if g in cred else "" for g in fm.index]
    fm["sl_expr_cohens_d"] = [round(cred.get(g, 0.0), 3) for g in fm.index]
    fm["sl_source"] = np.where(pd.Series(fm["sl_expr_cohens_d"].values, index=fm.index) > fm["sl_effect_v1"],
                               "expression_paralog", np.where(fm["sl_effect_v1"] > 0, "genotype", "none"))

# 3-tier novelty from tractability
ot_appr = tr["ot_sm_approved_drug"].astype(bool)
has_chem = tr["has_inhibitor_moa"].astype(bool) | (tr["chembl_n_potent"] >= 10)
tr["novelty_tier"] = np.where(ot_appr, "approved_sm_drug", np.where(has_chem, "has_sm_chemistry", "novel"))

# Expert-weighted composite (7 categories)
for c, cols in CAT_COLS.items():
    fm[f"cat_{c}"] = fm[cols].mean(axis=1)
fm["expert_composite"] = sum(fm[f"cat_{c}"] * W[c] for c in CATS)

# P2-SOFT base: all 283 ranked in one list, novelty tier as a visible tag
base = fm.sort_values("expert_composite", ascending=False).copy()
base["rank_all"] = range(1, len(base) + 1)
base["novelty_tier"] = tr["novelty_tier"].reindex(base.index)
gbm_onco_controls = {"PDGFRA", "MET", "FGFR1", "VRK1"}
base["is_gbm_oncogene_control"] = base.index.isin(gbm_onco_controls)
base.reset_index().to_csv(result(STEP, "nomination_base_p2soft.csv"), index=False)

# ---------------------------------------------------------------------------
# Build the top-100 evidence dossier (used by the live nomination prompt)
# ---------------------------------------------------------------------------
top100 = base.head(100).copy()


def _flags(g, r):
    f = []
    t = r["novelty_tier"]
    f.append({"novel": "NOVEL", "has_sm_chemistry": "has-chemistry", "approved_sm_drug": "APPROVED-DRUG"}[t])
    if not r.get("has_inhibitor_moa", False):
        f.append("no-inhibitor-MoA")
    if r.get("is_SL_partner_of_lost_TSG", False):
        f.append(f"SL:{r.get('lost_TSG_partner', '')}")
    elif isinstance(r.get("sl_best_context"), str) and r["sl_best_context"]:
        f.append(f"SL:{r['sl_best_context']}")
    if str(r.get("sc_compartment", "")) == "TME-dominated":
        f.append("TME-dominated")
    if r.get("cns_high_risk", False):
        f.append("CNS-safety-risk")
    if r.get("uniprot_kinase", False):
        f.append("kinase")
    elif r.get("uniprot_enzyme", False):
        f.append("enzyme")
    if r.get("is_common_essential", False):
        f.append("common-essential")
    return ";".join(x for x in f if x)


_lines = []
for i, (g, r) in enumerate(top100.iterrows(), 1):
    _lines.append(
        f"{i}. {g} | composite={r['expert_composite']:.3f} | "
        f"dep={r['cat_dependency']:.2f} selec={r['cat_selectivity']:.2f} SL={r['cat_synthetic_lethality']:.2f} "
        f"tumorSpec={r['cat_tumor_specificity']:.2f} tumorEnr={r['cat_tumor_enrichment']:.2f} "
        f"window={r['cat_therapeutic_window']:.2f} tract={r['cat_tractability']:.2f} | {_flags(g, r)}")
dossier = "\n".join(_lines)
genes = set(top100.index)

NOM_SYS = ("You are an expert in glioblastoma (GBM) small-molecule therapeutic-target discovery. "
           "You are reviewing a ranked list of candidate genes, each pre-scored on 7 evidence categories "
           "(dependency, selectivity vs pan-essential, synthetic-lethality with GBM lesions, tumor-cell specificity, "
           "tumor-vs-normal enrichment, therapeutic-window safety, small-molecule tractability). Composite is the "
           "expert-weighted sum. Your job: nominate the ~15 STRONGEST candidates as GBM small-molecule targets, "
           "integrating the scores with your biological judgment. Do not just copy the composite order — reason about "
           "mechanism, druggability, and directionality. A gene flagged 'no-inhibitor-MoA' or that functions as a tumor "
           "suppressor should be assessed for whether it needs INHIBITION vs REACTIVATION (small molecules cannot easily "
           "restore a lost function). Flag any candidate whose favorable score would be a mistake to act on.")
NOM_PROMPT = ("Candidate list (ranked by expert-weighted composite):\n\n" + dossier +
              "\n\nNominate the ~15 strongest GBM small-molecule TARGET candidates. For each: gene symbol, a one-sentence "
              "primary reason, the intended modality (inhibit / degrade / reactivate / unknown), and confidence "
              "(high/medium/low). Also list any genes you explicitly caution against despite a high score, with the reason.")
NOM_TOOL = {"name": "nominate_targets", "description": "Nominate the strongest GBM small-molecule target candidates.",
            "input_schema": {"type": "object", "properties": {
                "nominations": {"type": "array", "description": "~15 strongest candidates, best first",
                    "items": {"type": "object", "properties": {
                        "gene": {"type": "string"}, "reason": {"type": "string", "description": "one sentence"},
                        "modality": {"type": "string", "enum": ["inhibit", "degrade", "reactivate", "unknown"]},
                        "confidence": {"type": "string", "enum": ["high", "medium", "low"]}},
                        "required": ["gene", "reason", "modality", "confidence"]}},
                "cautions": {"type": "array", "description": "genes to caution against despite high score",
                    "items": {"type": "object", "properties": {"gene": {"type": "string"}, "concern": {"type": "string"}},
                        "required": ["gene", "concern"]}}},
                "required": ["nominations", "cautions"]}}


def nom_req():
    return {"prompt": NOM_PROMPT, "system": NOM_SYS, "tools": [NOM_TOOL],
            "tool_choice": {"type": "tool", "name": "nominate_targets"},
            "model": host.current_model(), "max_tokens": 2000}


# ---------------------------------------------------------------------------
# Nomination: cached runs, or live LLM ensemble
# ---------------------------------------------------------------------------
runs_path = result(STEP, "nomination_runs_raw.csv")
if use_cache(runs_path):
    print("Using cached nomination runs (no LLM calls).")
    nom = pd.read_csv(runs_path)
    n_err = 100 - nom["run"].nunique()
else:
    if not LIVE:
        raise FileNotFoundError(
            "Cached nomination_runs_raw.csv not found and GBM_LIVE is not set. Restore "
            "the committed file or set GBM_LIVE=1 (needs 'anthropic' + ANTHROPIC_API_KEY)."
        )
    print("LIVE: querying LLM nomination (100 runs)...")
    res = host.llm([nom_req() for _ in range(100)], max_concurrency=32)
    nom_rows, caution_rows, n_err = [], [], 0
    for run_i, r in enumerate(res):
        if "error" in r:
            n_err += 1
            continue
        inp = r["tool_use"]["input"]
        for rank, n in enumerate(inp.get("nominations", []), 1):
            g = n.get("gene", "")
            if g not in genes:
                continue
            nom_rows.append({"run": run_i, "nom_rank": rank, "gene": g, "modality": n.get("modality", ""),
                             "confidence": n.get("confidence", ""), "reason": n.get("reason", "")})
        for c in inp.get("cautions", []):
            caution_rows.append({"run": run_i, "gene": c.get("gene", ""), "concern": c.get("concern", "")})
    nom = pd.DataFrame(nom_rows)
    pd.DataFrame(caution_rows).to_csv(result(STEP, "nomination_cautions.csv"), index=False)
    nom.to_csv(runs_path, index=False)

# Nomination frequency table
freq = nom.groupby("gene").agg(n_runs=("run", "nunique"), mean_nom_rank=("nom_rank", "mean")).reset_index()
freq["nom_frequency"] = freq["n_runs"] / max(1, (100 - n_err))
mod = nom.groupby("gene")["modality"].agg(lambda s: s.value_counts().idxmax())
conf = nom.groupby("gene")["confidence"].agg(lambda s: s.value_counts().idxmax())
freq["dominant_modality"] = freq.gene.map(mod)
freq["dominant_confidence"] = freq.gene.map(conf)
freq = freq.sort_values(["nom_frequency", "mean_nom_rank"], ascending=[False, True])
freq.to_csv(result(STEP, "nomination_frequency.csv"), index=False)
freq_i = freq.set_index("gene")

# ---------------------------------------------------------------------------
# Final ranking table
# ---------------------------------------------------------------------------
final = base.copy()
final["nom_frequency"] = freq_i["nom_frequency"].reindex(final.index).fillna(0.0)
final["mean_nom_rank"] = freq_i["mean_nom_rank"].reindex(final.index)
final["dominant_modality"] = freq_i["dominant_modality"].reindex(final.index)
final["dominant_confidence"] = freq_i["dominant_confidence"].reindex(final.index)
final["llm_rank"] = final["nom_frequency"].rank(ascending=False, method="min").astype(int)

# Merge step11a gene-masked ablation columns when available.
masked_path = result("step11a_gene_masked_ablation", "gene_masked_ablation.csv")
if os.path.exists(masked_path):
    m = pd.read_csv(masked_path).set_index("gene")
    final["freq_named"] = m["freq_named"].reindex(final.index)
    final["freq_masked"] = m["freq_masked"].reindex(final.index)
    final["masked_delta"] = m["delta"].reindex(final.index)
    final["masked_effect_class"] = m["effect"].reindex(final.index).fillna("not-nominated")

gbm_onco = {"PDGFRA", "MET", "FGFR1"}


def _ctrl_class(g):
    if g in gbm_onco:
        return "GBM_oncogene_control"
    if base.loc[g, "novelty_tier"] == "approved_sm_drug":
        return "cytotoxic_or_panessential_drugged"
    return "candidate"


final["control_class"] = [_ctrl_class(g) for g in final.index]
final.reset_index().to_csv(result(STEP, "final_candidate_ranking.csv"), index=False)

print(f"\nFinal ranking: {len(final)} genes")
print("Top 10 by expert-weighted composite:")
for i, (g, r) in enumerate(final.head(10).iterrows(), 1):
    print(f"  {i:2d}. {g:10s} composite={r['expert_composite']:.3f} "
          f"nom_freq={r['nom_frequency']:.2f} tier={r['novelty_tier']}")
