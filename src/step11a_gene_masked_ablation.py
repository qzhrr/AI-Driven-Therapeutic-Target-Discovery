"""
step11a_gene_masked_ablation — Blinded (gene-masked) nomination ablation

The strongest test of whether the LLM nomination reflects the *evidence* or merely
recognizes *gene names*. Two arms of 100 nomination runs each over the same top-100
dossier:
  * named  — genes shown by symbol (identical to step10).
  * masked — gene symbols replaced by opaque codes (CX###); the numeric evidence
             (all 7 category scores + flags) is identical, only the identity is hidden.

A candidate whose nomination frequency holds up when masked is carried by its
evidence; one that collapses was riding name recognition. This is what separates
the genuinely novel, evidence-driven picks (VRK1, ELAVL1) from name-inflated ones.

Reproducibility note: the masking uses a fixed RNG seed (20260710), so the
code assignment is deterministic; the committed raw arms reproduce the published
effect classes exactly. A live re-run re-queries the LLM and is stochastic.

Outputs (results/step11a_gene_masked_ablation/):
  ablation_named_raw.csv, ablation_masked_raw.csv, gene_masked_ablation.csv

Run modes (see common.py):
  * cached (default): committed raw arms are used as-is — no LLM calls.
  * live (GBM_LIVE=1): re-query the LLM (needs 'anthropic' + ANTHROPIC_API_KEY).
"""
import os
import re
import sys
import json
import random

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import LIVE, use_cache, cache, result, ensure_dirs  # noqa: E402
from adapters import host  # noqa: E402

ensure_dirs()

STEP = "step11a_gene_masked_ablation"
N = 100
CATS = ["dependency", "selectivity", "synthetic_lethality", "tumor_specificity",
        "tumor_enrichment", "therapeutic_window", "tractability"]

# ---------------------------------------------------------------------------
# Inputs: dossier (committed cache) + base ranking (reconstructed from parquet)
# ---------------------------------------------------------------------------
d = json.load(open(cache("nom_dossier.json")))
dossier = d["dossier"]
genes = set(d["genes"])

# The step10 P2-SOFT base ranking, rebuilt from the committed feature matrix so
# this step does not depend on a transient intermediate.
fmx = pd.read_parquet(result("step08_feature_matrix", "feature_matrix.parquet"))
if "gene" in fmx.columns:
    fmx = fmx.set_index("gene")
fmx.index.name = "gene"
W = pd.read_csv(result("step09_weight_selection", "weights_aggregated.csv"))
wmap = dict(zip(W["category"], W["weight_normalized"]))
X = fmx[[f"cat_{c}" for c in CATS]].copy()
X.columns = CATS
fmx["expert_composite"] = X.values @ np.array([wmap[c] for c in CATS])

appr = fmx["ot_sm_approved_drug"] == True  # noqa: E712
sm_chem = (fmx["has_inhibitor_moa"] == True) | (fmx["chembl_n_potent"] >= 10)  # noqa: E712
tier = pd.Series("novel", index=fmx.index)
tier[sm_chem & ~appr] = "has_sm_chemistry"
tier[appr] = "approved_sm_drug"
fmx["novelty_tier"] = tier
fmx = fmx.sort_values("expert_composite", ascending=False)
fmx["rank_all"] = np.arange(1, len(fmx) + 1)
fmx["is_gbm_oncogene_control"] = fmx.index.isin({"PDGFRA", "MET", "FGFR1", "VRK1"})
final = fmx.reset_index()[["gene", "novelty_tier", "rank_all",
                           "cat_selectivity", "cat_dependency", "cat_tractability",
                           "is_gbm_oncogene_control"]]

# ---------------------------------------------------------------------------
# Build the masked dossier (deterministic code assignment)
# ---------------------------------------------------------------------------
lines = [ln for ln in dossier.split("\n") if ln.strip()]
pat = re.compile(r"^(\d+)\.\s+([^|]+?)\s+\|\s+(.*)$")
parsed = []
for ln in lines:
    m = pat.match(ln)
    assert m, f"unparsed: {ln[:80]}"
    parsed.append((int(m.group(1)), m.group(2).strip(), m.group(3)))

rng = random.Random(20260710)
codes = [f"CX{n:03d}" for n in rng.sample(range(100, 1000), len(parsed))]
code2gene, masked_lines = {}, []
for (rk, g, rest), code in zip(parsed, codes):
    code2gene[code] = g
    assert g not in rest, f"gene leaks into flags: {g}"
    masked_lines.append(f"{rk}. {code} | {rest}")
masked_dossier = "\n".join(masked_lines)
codes_set = set(code2gene.keys())

# ---------------------------------------------------------------------------
# Prompts / tool (named and masked framings)
# ---------------------------------------------------------------------------
NOM_SYS = ("You are an expert in glioblastoma (GBM) small-molecule therapeutic-target discovery. "
           "You are reviewing a ranked list of candidate genes, each pre-scored on 7 evidence categories "
           "(dependency, selectivity vs pan-essential, synthetic-lethality with GBM lesions, tumor-cell specificity, "
           "tumor-vs-normal enrichment, therapeutic-window safety, small-molecule tractability). Composite is the "
           "expert-weighted sum. Your job: nominate the ~15 STRONGEST candidates as GBM small-molecule targets, "
           "integrating the scores with your biological judgment. Do not just copy the composite order — reason about "
           "mechanism, druggability, and directionality. A gene flagged 'no-inhibitor-MoA' or that functions as a tumor "
           "suppressor should be assessed for whether it needs INHIBITION vs REACTIVATION (small molecules cannot easily "
           "restore a lost function). Flag any candidate whose favorable score would be a mistake to act on.")
NOM_SYS_M = NOM_SYS.replace(
    "reviewing a ranked list of candidate genes",
    "reviewing a ranked list of ANONYMIZED candidates (identified only by opaque codes like CX123; gene names are withheld)")

NOM_TOOL = {"name": "nominate_targets", "description": "Nominate the strongest GBM small-molecule target candidates.",
            "input_schema": {"type": "object", "properties": {
                "nominations": {"type": "array", "description": "~15 strongest candidates, best first",
                    "items": {"type": "object", "properties": {
                        "gene": {"type": "string"}, "reason": {"type": "string", "description": "one sentence"},
                        "modality": {"type": "string", "enum": ["inhibit", "degrade", "reactivate", "unknown"]},
                        "confidence": {"type": "string", "enum": ["high", "medium", "low"]}},
                        "required": ["gene", "reason", "modality", "confidence"]}},
                "cautions": {"type": "array", "description": "candidates to caution against despite high score",
                    "items": {"type": "object", "properties": {"gene": {"type": "string"}, "concern": {"type": "string"}},
                        "required": ["gene", "concern"]}}},
                "required": ["nominations", "cautions"]}}


def prompt_for(doss, ident_word):
    return ("Candidate list (ranked by expert-weighted composite):\n\n" + doss +
            f"\n\nNominate the ~15 strongest GBM small-molecule TARGET candidates. For each: {ident_word}, a one-sentence "
            "primary reason, the intended modality (inhibit / degrade / reactivate / unknown), and confidence "
            "(high/medium/low). Also list any candidates you explicitly caution against despite a high score, with the reason.")


def req(prompt, sysp):
    return {"prompt": prompt, "system": sysp, "tools": [NOM_TOOL],
            "tool_choice": {"type": "tool", "name": "nominate_targets"},
            "model": host.current_model(), "max_tokens": 2000}


def run_arm(prompt, sysp, valid, arm):
    res = host.llm([req(prompt, sysp) for _ in range(N)], max_concurrency=32)
    rows, err = [], 0
    for i, r in enumerate(res):
        if "error" in r:
            err += 1
            continue
        inp = r["tool_use"]["input"]
        for rk, n in enumerate(inp.get("nominations", []), 1):
            gid = n.get("gene", "").strip()
            if gid not in valid:
                continue
            rows.append({"arm": arm, "run": i, "nom_rank": rk, "ident": gid,
                         "modality": n.get("modality", ""), "confidence": n.get("confidence", ""),
                         "reason": n.get("reason", "")})
    return pd.DataFrame(rows), err


# ---------------------------------------------------------------------------
# Obtain both arms (cached, or live)
# ---------------------------------------------------------------------------
named_path = result(STEP, "ablation_named_raw.csv")
masked_path = result(STEP, "ablation_masked_raw.csv")

if use_cache(named_path) and use_cache(masked_path):
    print("Using cached ablation arms (no LLM calls).")
    named_df = pd.read_csv(named_path)
    masked_df = pd.read_csv(masked_path)
else:
    if not LIVE:
        raise FileNotFoundError(
            "Cached ablation arms not found and GBM_LIVE is not set. Restore the "
            "committed ablation_*_raw.csv or set GBM_LIVE=1 (needs 'anthropic' + ANTHROPIC_API_KEY)."
        )
    print(f"LIVE: querying LLM ablation ({N} runs x 2 arms)...")
    named_df, _ = run_arm(prompt_for(dossier, "gene symbol"), NOM_SYS, genes, "named")
    masked_df, _ = run_arm(prompt_for(masked_dossier, "the candidate code (e.g. CX123) in the 'gene' field"),
                           NOM_SYS_M, codes_set, "masked")
    masked_df["gene"] = masked_df["ident"].map(code2gene)
    named_df["gene"] = named_df["ident"]
    named_df.to_csv(named_path, index=False)
    masked_df.to_csv(masked_path, index=False)

# ---------------------------------------------------------------------------
# Named vs masked frequency + effect classification
# ---------------------------------------------------------------------------
n_named, n_masked = named_df["run"].nunique(), masked_df["run"].nunique()


def freq_table(df, nruns):
    g = df.groupby("gene").agg(n_runs=("run", "nunique"), mean_rank=("nom_rank", "mean")).reset_index()
    g["freq"] = g["n_runs"] / nruns
    return g.set_index("gene")


fn = freq_table(named_df, n_named)
fmk = freq_table(masked_df, n_masked)
comp = pd.DataFrame(index=sorted(set(fn.index) | set(fmk.index)))
comp["freq_named"] = fn["freq"].reindex(comp.index).fillna(0)
comp["freq_masked"] = fmk["freq"].reindex(comp.index).fillna(0)
comp["delta"] = comp["freq_masked"] - comp["freq_named"]

final_idx = final.set_index("gene")
for c in ["novelty_tier", "rank_all", "cat_selectivity", "cat_dependency", "cat_tractability", "is_gbm_oncogene_control"]:
    comp[c] = final_idx[c].reindex(comp.index)


def cls(r):
    if r["freq_named"] >= 0.5 and r["freq_masked"] >= 0.5 * r["freq_named"]:
        return "evidence-carried (robust)"
    if r["freq_named"] >= 0.5 and r["freq_masked"] < 0.5 * r["freq_named"]:
        return "NAME-INFLATED (drops when masked)"
    if r["freq_named"] < 0.3 and r["freq_masked"] >= 0.5:
        return "NAME-SUPPRESSED (rises when masked)"
    return "other/low"


act = comp[(comp["freq_named"] >= 0.10) | (comp["freq_masked"] >= 0.10)].copy()
act["effect"] = act.apply(cls, axis=1)

out = comp.copy()
out.index.name = "gene"
out = out.reset_index()
out["effect"] = out.apply(lambda r: (act.loc[r["gene"], "effect"] if r["gene"] in act.index else "not-nominated"), axis=1)
out = out[["gene", "freq_named", "freq_masked", "delta", "effect", "novelty_tier", "rank_all",
           "cat_selectivity", "cat_dependency", "cat_tractability", "is_gbm_oncogene_control"]]
out = out.sort_values("freq_named", ascending=False)
out.to_csv(result(STEP, "gene_masked_ablation.csv"), index=False)

print(f"\nGene-masked ablation: {len(out)} genes")
print("Most name-inflated (largest drop when masked):")
drops = out[out["effect"].str.contains("NAME-INFLATED", na=False)].nlargest(5, "freq_named")
for _, r in drops.iterrows():
    print(f"  {r['gene']:10s} named={r['freq_named']:.2f} masked={r['freq_masked']:.2f}")
print("Robust novel picks (evidence-carried, novel tier):")
rob = out[(out["effect"].str.contains("robust", na=False)) & (out["novelty_tier"] == "novel")].nlargest(5, "freq_masked")
for _, r in rob.iterrows():
    print(f"  {r['gene']:10s} named={r['freq_named']:.2f} masked={r['freq_masked']:.2f}")
