"""
step09_weight_selection — LLM weight selection (single-model ensemble)

Distributes 100 points across the 7 evidence categories using an LLM, under two
framings run 100 times each:
  * expert     — GBM small-molecule target-discovery system prompt + domain category
                 definitions.
  * non_expert — generic feature framing with all domain cues stripped (control for
                 how much the domain framing, not the data, drives the weights).

The final weight vector is the per-category expert median, normalized to sum 1.
Weights are assigned BLIND to candidate scores (domain reasoning only), which
firewalls the weighting from the ranking it later produces.

Outputs (results/step09_weight_selection/):
  weights_expert_raw.csv, weights_nonexpert_raw.csv, weights_aggregated.csv,
  weights_comparison.csv   and results/.../weight_vector.json (consumed by step10).

Run modes (see common.py):
  * cached (default): committed outputs are used as-is — no LLM calls. The
    published weight vector is reproduced exactly.
  * live (GBM_LIVE=1): re-query the LLM (needs the 'anthropic' package and
    ANTHROPIC_API_KEY). Being stochastic, a live re-run draws new samples and
    will not reproduce the committed values bit-for-bit.
"""
import os
import sys
import json

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import LIVE, use_cache, result, ensure_dirs  # noqa: E402
from adapters import host  # noqa: E402

ensure_dirs()

STEP = "step09_weight_selection"
N = 100  # runs per condition

CATS = ["dependency", "selectivity", "synthetic_lethality", "tumor_specificity",
        "tumor_enrichment", "therapeutic_window", "tractability"]

EXPERT_DEFS = {
    "dependency": "How selectively glioma cells depend on this gene for survival (CRISPR knockout reduces fitness in glioma lines more than other lineages) — direct genetic evidence the tumor needs it.",
    "selectivity": "Penalty for pan-essential genes required by virtually all cell types; high weight favors glioma-selective dependencies over generic housekeeping requirements (which are toxic to inhibit).",
    "synthetic_lethality": "Whether the dependency is synthetic-lethal with a recurrent GBM lesion (CDKN2A/MTAP deletion, PTEN/RB1/NF1 loss, EGFR/PDGFRA amplification) — a genotype-stratified vulnerability.",
    "tumor_specificity": "Whether the gene is expressed specifically in malignant tumor cells rather than the non-malignant microenvironment (immune/vascular/normal brain), from single-cell data.",
    "tumor_enrichment": "Whether the gene is over-expressed in GBM tumor tissue relative to normal brain across a patient cohort.",
    "therapeutic_window": "Safety: how low the gene's expression is in normal brain and vital organs; high weight favors a wide therapeutic window with less on-target toxicity.",
    "tractability": "How amenable the target is to small-molecule drugging — druggable pocket, ligand-bound structure, druggable family, existing chemical matter.",
}
NONEXPERT_DEFS = {
    "dependency": "A score measuring how strongly one group of samples shows a signal relative to other groups.",
    "selectivity": "A score penalizing features that are uniformly high across all groups.",
    "synthetic_lethality": "A score measuring association between a feature and a categorical subgroup label.",
    "tumor_specificity": "A score measuring whether a value is concentrated in one category versus others.",
    "tumor_enrichment": "A score comparing values between two sample groups.",
    "therapeutic_window": "A score measuring how low a value is across a set of reference categories.",
    "tractability": "A score measuring a structural/annotation property of the entity.",
}
EXPERT_SYS = ("You are an expert in glioblastoma (GBM) small-molecule therapeutic-target discovery. "
             "You are setting importance weights for categories of evidence used to prioritize candidate genes "
             "as druggable small-molecule inhibition targets in GBM.")
NONEXPERT_SYS = ("You are assigning importance weights to categories of numerical features in a dataset. "
                 "Base your allocation only on the category names and definitions given.")
TOOL = {"name": "assign_weights",
        "description": "Assign integer importance weights (points, summing to ~100) to the 7 categories.",
        "input_schema": {"type": "object", "properties": {
            **{c: {"type": "integer", "description": "points 0-100"} for c in CATS},
            "rationale": {"type": "string", "description": "one concise sentence justifying the allocation"}},
            "required": CATS + ["rationale"]}}


def build_prompt(defs):
    lines = [f"- {c}: {defs[c]}" for c in CATS]
    return ("Distribute exactly 100 points across these 7 categories to reflect how important each should be "
            "in the final weighted score. More important = more points. Use integers; they should sum to 100.\n\n"
            + "\n".join(lines))


def make_req(cond):
    defs = EXPERT_DEFS if cond == "expert" else NONEXPERT_DEFS
    sysp = EXPERT_SYS if cond == "expert" else NONEXPERT_SYS
    return {"prompt": build_prompt(defs), "system": sysp, "tools": [TOOL],
            "tool_choice": {"type": "tool", "name": "assign_weights"},
            "model": host.current_model(), "max_tokens": 500}


def run_condition(cond):
    """Return a normalized (points scaled to sum 100) per-run DataFrame."""
    res = host.llm([make_req(cond) for _ in range(N)], max_concurrency=32)
    rows, err = [], 0
    for i, r in enumerate(res):
        if "error" in r:
            err += 1
            continue
        w = r["tool_use"]["input"]
        row = {"run": i, "condition": cond}
        row.update({c: w[c] for c in CATS})
        row["rationale"] = w.get("rationale", "")
        row["sum"] = sum(w[c] for c in CATS)
        rows.append(row)
    df = pd.DataFrame(rows)
    print(f"  {cond}: {len(df)} ok, {err} errors")
    for c in CATS:
        df[c] = df[c] / df["sum"] * 100
    return df


# ---------------------------------------------------------------------------
# 1. Obtain per-run weights (from cache, or by querying the LLM live)
# ---------------------------------------------------------------------------
exp_path = result(STEP, "weights_expert_raw.csv")
ne_path = result(STEP, "weights_nonexpert_raw.csv")

if use_cache(exp_path) and use_cache(ne_path):
    print("Using cached raw weight runs (no LLM calls).")
    exp_df = pd.read_csv(exp_path)
    ne_df = pd.read_csv(ne_path)
    # Committed files store already-normalized weights; guard re-normalization.
    for df in (exp_df, ne_df):
        if "sum" in df.columns and not np.allclose(df[CATS].sum(axis=1), 100, atol=1):
            for c in CATS:
                df[c] = df[c] / df["sum"] * 100
else:
    if not LIVE:
        raise FileNotFoundError(
            "Cached weight runs not found and GBM_LIVE is not set. Either restore "
            "the committed results/step09_weight_selection/weights_*_raw.csv or set "
            "GBM_LIVE=1 (needs the 'anthropic' package + ANTHROPIC_API_KEY) to recompute."
        )
    print(f"LIVE: querying LLM ({N} runs x 2 conditions)...")
    exp_df = run_condition("expert")
    ne_df = run_condition("non_expert")
    exp_df.to_csv(exp_path, index=False)
    ne_df.to_csv(ne_path, index=False)

# ---------------------------------------------------------------------------
# 2. Aggregate expert runs -> final weight vector
# ---------------------------------------------------------------------------
med = exp_df[CATS].median()
wvec = (med / med.sum()).round(4)
agg = pd.DataFrame({
    "category": CATS,
    "expert_median_points": med.values,
    "weight_normalized": wvec.values,
    "expert_mean": exp_df[CATS].mean().values,
    "expert_std": exp_df[CATS].std().values,
    "expert_cv": (exp_df[CATS].std() / exp_df[CATS].mean()).values,
}).sort_values("weight_normalized", ascending=False)
agg.to_csv(result(STEP, "weights_aggregated.csv"), index=False)
json.dump({c: float(w) for c, w in zip(CATS, wvec.values)},
          open(result(STEP, "weight_vector.json"), "w"))

print("\nFINAL WEIGHT VECTOR (expert median, normalized to 1.0):")
for _, r in agg.iterrows():
    print(f"  {r['category']:20s} {r['weight_normalized']:.3f}  (CV {r['expert_cv']:.3f})")
print("  sum:", round(float(wvec.sum()), 4))
print("  mean CV across categories:",
      round(float((exp_df[CATS].std() / exp_df[CATS].mean()).mean()), 3))

# ---------------------------------------------------------------------------
# 3. Expert vs non-expert comparison
# ---------------------------------------------------------------------------
comp = []
for c in CATS:
    e, n = exp_df[c], ne_df[c]
    comp.append({"category": c,
                 "expert_median": e.median(), "expert_mean": e.mean(), "expert_std": e.std(),
                 "nonexpert_median": n.median(), "nonexpert_mean": n.mean(), "nonexpert_std": n.std(),
                 "delta": e.median() - n.median(), "mw_p": mannwhitneyu(e, n).pvalue})
pd.DataFrame(comp).to_csv(result(STEP, "weights_comparison.csv"), index=False)
print("\nWrote weights_aggregated.csv, weights_comparison.csv, weight_vector.json")
