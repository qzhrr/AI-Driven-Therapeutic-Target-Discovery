"""
step04_tumor_vs_tme — Tumor-vs-TME deconvolution

Pseudobulk GBmap single-cell; malignant-vs-microenvironment specificity per candidate.

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

# Tumor-vs-TME deconvolution (GBmap single-cell)


import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

# Load data
A = ad.read_h5ad(raw("core_gbmap.h5ad"), backed="r")

candidate_genes_df = pd.read_csv(result("step01_dependency", "glioma_candidates.csv"), index_col="gene")
candidate_genes = candidate_genes_df.index.tolist()

# Map symbol -> ensembl via feature_name
var = A.var.copy()
sym2ens = {}
for ens, row in zip(var.index, var["feature_name"]):
    sym2ens.setdefault(str(row), ens)
found = {g: sym2ens[g] for g in candidate_genes if g in sym2ens}
ens_ids = list(found.values())
gene_syms = list(found.keys())

# Column indices for our genes (in var order)
var_pos = {e: i for i, e in enumerate(A.var.index)}
col_idx = np.array([var_pos[e] for e in ens_ids])

cell_type = A.obs["cell_type"].astype(str).values
types = pd.unique(cell_type)

# Accumulators: sum of expression and count of expressing cells, per (type x gene)
n_g = len(col_idx)
sum_expr = {t: np.zeros(n_g) for t in types}
cnt_expr = {t: np.zeros(n_g) for t in types}
n_cells = {t: 0 for t in types}

CHUNK = 20000
N = A.shape[0]
for start in range(0, N, CHUNK):
    end = min(start + CHUNK, N)
    Xc = A.X[start:end]
    Xc = Xc[:, col_idx]
    if sparse.issparse(Xc):
        Xc = Xc.toarray()
    Xc = np.asarray(Xc)
    ct_chunk = cell_type[start:end]
    for t in np.unique(ct_chunk):
        m = ct_chunk == t
        sub = Xc[m]
        sum_expr[t] += sub.sum(0)
        cnt_expr[t] += (sub > 0).sum(0)
        n_cells[t] += int(m.sum())

mean_expr = pd.DataFrame({t: sum_expr[t] / max(n_cells[t], 1) for t in types}, index=gene_syms)
pct_expr = pd.DataFrame({t: cnt_expr[t] / max(n_cells[t], 1) for t in types}, index=gene_syms)

MALIG = "malignant cell"
tme_types = [t for t in mean_expr.columns if t != MALIG]
robust_tme = [t for t in tme_types if n_cells[t] >= 500]

malig = mean_expr[MALIG]
tme_max = mean_expr[robust_tme].max(axis=1)
tme_max_type = mean_expr[robust_tme].idxmax(axis=1)
tme_mean = mean_expr[robust_tme].mean(axis=1)

sc = pd.DataFrame(index=mean_expr.index)
sc["sc_malignant_mean"] = malig
sc["sc_tme_max_mean"] = tme_max
sc["sc_tme_max_type"] = tme_max_type
sc["sc_tme_mean"] = tme_mean
sc["sc_pct_malignant_expressing"] = pct_expr[MALIG]
sc["sc_malignant_specificity"] = malig - tme_max
sc["sc_tme_dominated_flag"] = (tme_max > malig + 0.25)

# Add collapsed malignant-vs-all-TME contrast using level1
lvl1 = A.obs["annotation_level_1"].astype(str).values
sum_neo = np.zeros(len(col_idx))
sum_non = np.zeros(len(col_idx))
cnt_neo = np.zeros(len(col_idx))
cnt_non = np.zeros(len(col_idx))
nneo = 0
nnon = 0

CHUNK = 40000
for s in range(0, A.shape[0], CHUNK):
    e = min(s + CHUNK, A.shape[0])
    Xc = A.X[s:e][:, col_idx]
    Xc = Xc.toarray() if sparse.issparse(Xc) else np.asarray(Xc)
    l = lvl1[s:e]
    mn = l == "Neoplastic"
    mo = l == "Non-neoplastic"
    sum_neo += Xc[mn].sum(0)
    cnt_neo += (Xc[mn] > 0).sum(0)
    nneo += int(mn.sum())
    sum_non += Xc[mo].sum(0)
    cnt_non += (Xc[mo] > 0).sum(0)
    nnon += int(mo.sum())

sc["sc_neoplastic_mean"] = sum_neo / nneo
sc["sc_nonneoplastic_mean"] = sum_non / nnon
sc["sc_malignant_vs_tme_collapsed"] = sc["sc_neoplastic_mean"] - sc["sc_nonneoplastic_mean"]

# Reindex to all 283 candidates (unmapped -> NaN)
full = pd.DataFrame(index=candidate_genes)
sc_full = full.join(sc)
sc_full["sc_measured"] = sc_full["sc_malignant_mean"].notna()
sc_full.to_csv(result("step04_tumor_vs_tme", "singlecell_compartment_expression.csv"))


import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

# Load the GBmap atlas
A = ad.read_h5ad(raw("core_gbmap.h5ad"), backed="r")

# Load candidate genes
candidate_genes = pd.read_csv(result("step01_dependency", "glioma_candidates.csv"), index_col="gene").index.tolist()

# Map symbol -> ensembl via feature_name
var = A.var.copy()
sym2ens = {}
for ens, row in zip(var.index, var["feature_name"]):
    sym2ens.setdefault(str(row), ens)
found = {g: sym2ens[g] for g in candidate_genes if g in sym2ens}
ens_ids = list(found.values())
gene_syms = list(found.keys())

# Column indices for our genes (in var order)
var_pos = {e: i for i, e in enumerate(A.var.index)}
col_idx = np.array([var_pos[e] for e in ens_ids])

cell_type = A.obs["cell_type"].astype(str).values
types = pd.unique(cell_type)

# Accumulators: sum of expression and count of expressing cells, per (type x gene)
n_g = len(col_idx)
sum_expr = {t: np.zeros(n_g) for t in types}
cnt_expr = {t: np.zeros(n_g) for t in types}
n_cells = {t: 0 for t in types}

CHUNK = 20000
N = A.shape[0]
for start in range(0, N, CHUNK):
    end = min(start + CHUNK, N)
    Xc = A.X[start:end]
    Xc = Xc[:, col_idx]
    if sparse.issparse(Xc):
        Xc = Xc.toarray()
    Xc = np.asarray(Xc)
    ct_chunk = cell_type[start:end]
    for t in np.unique(ct_chunk):
        m = ct_chunk == t
        sub = Xc[m]
        sum_expr[t] += sub.sum(0)
        cnt_expr[t] += (sub > 0).sum(0)
        n_cells[t] += int(m.sum())
    if start % 100000 == 0:
        print(f"  {end}/{N}")

mean_expr = pd.DataFrame({t: sum_expr[t] / max(n_cells[t], 1) for t in types}, index=gene_syms)
pct_expr = pd.DataFrame({t: cnt_expr[t] / max(n_cells[t], 1) for t in types}, index=gene_syms)

mean_expr.to_csv(result("step04_tumor_vs_tme", "gbmap_pseudobulk_mean.csv"))
print("saved gbmap_pseudobulk_mean.csv:", mean_expr.shape)
