"""
step05_therapeutic_window — Therapeutic window / safety

Graded percentile safety over GTEx + brain subregions + Lake normal-brain neurons.

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

# Therapeutic window / safety (GTEx + Lake + single-cell)


import json, numpy as np, pandas as pd

main = json.load(open(cache("gtex_median.json")))
alias = json.load(open(cache("gtex_median_alias.json")))
sym_map = {**main["sym_by_gencode"], **alias["sym_by_gencode"]}
medians = main["rows"]["medians"] + alias["rows"]

recs = []
for m in medians:
    sym = sym_map.get(m["gencodeId"], m["geneSymbol"])
    recs.append((sym, m["tissueSiteDetailId"], m["median"]))
gtex = pd.DataFrame(recs, columns=["gene", "tissue", "tpm"]).pivot_table(index="gene", columns="tissue", values="tpm", aggfunc="max")

brain_tissues = [t for t in gtex.columns if t.startswith("Brain")]
vital = ["Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Liver", "Lung", "Kidney_Cortex", "Kidney_Medulla"]
vital = [v for v in vital if v in gtex.columns]
nonbrain = [t for t in gtex.columns if not t.startswith("Brain")]

feat = pd.DataFrame(index=gtex.index)
feat["gtex_max_tpm"] = gtex.max(axis=1)
feat["gtex_median_tpm"] = gtex.median(axis=1)
feat["gtex_max_tissue"] = gtex.idxmax(axis=1)
feat["gtex_vital_organ_max"] = gtex[vital].max(axis=1)
feat["gtex_brain_max"] = gtex[brain_tissues].max(axis=1)
feat["gtex_brain_mean"] = gtex[brain_tissues].mean(axis=1)
feat["gtex_nonbrain_max"] = gtex[nonbrain].max(axis=1)
feat.to_csv(result("step05_therapeutic_window", "gtex_safety_features.csv"))


import gzip, numpy as np, pandas as pd, re

cand = pd.read_csv(result("step01_dependency", "glioma_candidates.csv"))["gene"].tolist()
cand_set = set(cand)

def celltype_class(colname):
    p = colname.strip('"').split("_")[0]
    if re.match(r"^Ex\d", p): return "Excitatory_neuron"
    if re.match(r"^In\d", p): return "Inhibitory_neuron"
    if p in ("Gran", "Purk1", "Purk2"): return "Cerebellar_neuron"
    return {"Ast": "Astrocyte", "Oli": "Oligodendrocyte", "OPC": "OPC", "Mic": "Microglia",
            "End": "Endothelial", "Per": "Pericyte"}.get(p, p)

files = {
    "FrontalCortex": raw("GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz"),
    "VisualCortex": raw("GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz"),
    "CerebellarHem": raw("GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz"),
}

region_ct_sum = {}
for region, path in files.items():
    with gzip.open(path, "rt") as f:
        header = f.readline()
        ncells = len(header.rstrip("\n").split("\t"))
        cells = [c.strip('"') for c in header.rstrip("\n").split("\t")]
    classes = [celltype_class(c) for c in cells]
    with gzip.open(path, "rt") as f:
        f.readline()
        mat_rows = []; genes_kept = []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            g = parts[0].strip('"')
            if g in cand_set:
                mat_rows.append(np.array(parts[1:], dtype=np.float32))
                genes_kept.append(g)
    M = np.vstack(mat_rows)
    print(f"{region}: {M.shape[0]} candidate genes x {len(cells)} cells")
    region_ct_sum[region] = (M, genes_kept, np.array(classes))
print("parsed all 3 regions")

libsize = {}
for region, path in files.items():
    with gzip.open(path, "rt") as f:
        header = f.readline()
        ncells = len(header.rstrip("\n").split("\t"))
        tot = np.zeros(ncells, dtype=np.float64)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            tot += np.array(parts[1:], dtype=np.float64)
    libsize[region] = tot
    print(f"{region}: median libsize {np.median(tot):.0f}, cells {ncells}")

def agg(region):
    M, genes_kept, classes = region_ct_sum[region]
    ls = libsize[region]
    cpm = M / ls[None, :] * 1e4
    logcpm = np.log1p(cpm)
    df = pd.DataFrame(logcpm, index=genes_kept, columns=classes)
    return df.T.groupby(level=0).mean().T

per_region = {r: agg(r) for r in files}
for r, df in per_region.items():
    print(f"\n{r}: {df.shape[0]} genes x {df.shape[1]} celltypes -> {list(df.columns)}")

neuron_cols = {
    "FrontalCortex": ["Excitatory_neuron", "Inhibitory_neuron"],
    "VisualCortex": ["Excitatory_neuron", "Inhibitory_neuron"],
    "CerebellarHem": ["Cerebellar_neuron"],
}

rows = []
for g in cand:
    neuron_vals = []; any_vals = []
    for r, df in per_region.items():
        if g in df.index:
            ncs = [c for c in neuron_cols[r] if c in df.columns]
            neuron_vals += [df.loc[g, c] for c in ncs]
            any_vals += [df.loc[g, c] for c in df.columns]
    rows.append({"gene": g,
                 "lake_neuron_max": max(neuron_vals) if neuron_vals else np.nan,
                 "lake_neuron_mean": np.mean(neuron_vals) if neuron_vals else np.nan,
                 "lake_brain_max": max(any_vals) if any_vals else np.nan,
                 "lake_measured": len(any_vals) > 0})
lake = pd.DataFrame(rows).set_index("gene")
print("\nLake features built:", lake.shape, "| measured:", int(lake.lake_measured.sum()))
print("\nHighest normal-neuron expression (top CNS-toxicity risk):")
print(lake.sort_values("lake_neuron_max", ascending=False).head(10)[["lake_neuron_max", "lake_brain_max"]].round(2).to_string())
lake.to_csv(result("step05_therapeutic_window", "lake2018_brain_features.csv"))


import pandas as pd
import numpy as np

# Apply figure style (global state)
apply_figure_style()

# Load inputs
cand = pd.read_csv(result("step01_dependency", "glioma_candidates.csv")).set_index("gene")
sl = pd.read_csv(result("step03_synthetic_lethality", "sl_candidate_features.csv"), index_col=0)
sc = pd.read_csv(result("step04_tumor_vs_tme", "singlecell_compartment_expression.csv"), index_col=0)
gbmap_pb = pd.read_csv(result("step04_tumor_vs_tme", "gbmap_pseudobulk_mean.csv"), index_col=0)

# Load GTEx and Lake features (produced earlier in the session, saved locally)
gtex = pd.read_csv(result("step05_therapeutic_window", "gtex_safety_features.csv"), index_col=0)
lake = pd.read_csv(result("step05_therapeutic_window", "lake2018_brain_features.csv"), index_col=0)

# GBmap normal-brain compartments
gbmap_normal_brain_cols = [c for c in ["oligodendrocyte", "astrocyte", "radial glial cell", "neuron",
                                        "oligodendrocyte precursor cell"] if c in gbmap_pb.columns]
gbmap_norm_brain = gbmap_pb[gbmap_normal_brain_cols].max(axis=1)

tw = pd.DataFrame(index=cand.index)
# normal tissue (body-wide)
tw["gtex_max_tpm"] = gtex["gtex_max_tpm"]
tw["gtex_vital_organ_max"] = gtex["gtex_vital_organ_max"]
tw["gtex_max_tissue"] = gtex["gtex_max_tissue"]
# normal brain (3 sources)
tw["gtex_brain_max"] = gtex["gtex_brain_max"]
tw["lake_neuron_max"] = lake["lake_neuron_max"]
tw["lake_brain_max"] = lake["lake_brain_max"]
tw["gbmap_normal_brain_max"] = gbmap_norm_brain
# tumor / malignant expression (from Step 4)
tw["sc_malignant_mean"] = sc["sc_malignant_mean"]
tw["sc_tme_dominated_flag"] = sc["sc_tme_dominated_flag"]
tw["sc_malignant_specificity"] = sc["sc_malignant_specificity"]

# CNS-safety flag
lake_neuron_hi = lake["lake_neuron_max"].quantile(0.75)
tw["cns_safety_flag"] = ((tw["gtex_brain_max"] > 10) | (tw["lake_neuron_max"] > lake_neuron_hi)).fillna(False)
# vital-organ safety flag
tw["vital_organ_flag"] = (tw["gtex_vital_organ_max"] > 25).fillna(False)

# Graded safety scores
def safety_score(series, log=True):
    s = series.copy()
    if log:
        s = np.log1p(s)
    r = s.rank(pct=True)
    return 1.0 - r

tw["safety_vital_organ"] = safety_score(tw["gtex_vital_organ_max"])
tw["safety_brain_bulk"] = safety_score(tw["gtex_brain_max"])
tw["safety_brain_neuron"] = safety_score(tw["lake_neuron_max"].fillna(tw["lake_neuron_max"].median()))
tw["safety_body_wide"] = safety_score(tw["gtex_max_tpm"])
tw["safety_composite"] = tw[["safety_vital_organ", "safety_brain_bulk", "safety_brain_neuron"]].mean(axis=1)

# Hard flags for extreme cases (top-decile risk)
tw["cns_high_risk"] = (tw["gtex_brain_max"] >= tw["gtex_brain_max"].quantile(0.90)) | \
                      (tw["lake_neuron_max"] >= tw["lake_neuron_max"].quantile(0.90))
tw["vital_high_risk"] = tw["gtex_vital_organ_max"] >= tw["gtex_vital_organ_max"].quantile(0.90)

# Reviewer item 1: min_subgroup_n fragility field
n_glioma = 72
tw["glioma_dep_fraction"] = cand["glioma_dep_fraction"]
tw["archetype"] = cand["archetype"]
tw["min_subgroup_n"] = (cand["glioma_dep_fraction"] * n_glioma).round().astype(int)
if "sl_n_contexts" in sl.columns:
    tw["sl_n_contexts"] = sl["sl_n_contexts"].reindex(tw.index).fillna(0).astype(int)

# Reviewer item 2: sc_specificity_tier
def spec_tier(r):
    if r.get("sc_tme_dominated_flag") == True:
        return "TME-dominated"
    spec = r.get("sc_malignant_specificity")
    if pd.isna(spec):
        return "not-measured"
    if spec >= 0.15:
        return "confirmed-malignant-enriched"
    return "not-TME-dominated-weak"

tw["sc_specificity_tier"] = sc.reindex(tw.index).apply(spec_tier, axis=1)
tw["sc_compartment"] = np.where(sc["sc_tme_dominated_flag"].reindex(tw.index) == True,
                                "TME-dominated", "tumor-cell-expressed")

# Reorder columns and save
tw_out = tw.copy()
tw_out.index.name = "gene"
col_order = ["archetype", "glioma_dep_fraction", "min_subgroup_n", "sl_n_contexts",
             "sc_compartment", "sc_specificity_tier", "sc_malignant_mean", "sc_malignant_specificity", "sc_tme_dominated_flag",
             "gtex_max_tpm", "gtex_max_tissue", "gtex_vital_organ_max", "gtex_brain_max",
             "lake_neuron_max", "lake_brain_max", "gbmap_normal_brain_max",
             "safety_vital_organ", "safety_brain_bulk", "safety_brain_neuron", "safety_body_wide", "safety_composite",
             "cns_high_risk", "vital_high_risk"]
col_order = [c for c in col_order if c in tw_out.columns]
tw_out = tw_out[col_order]
tw_out.to_csv(result("step05_therapeutic_window", "therapeutic_window_features.csv"))
print("saved therapeutic_window_features.csv:", tw_out.shape)
