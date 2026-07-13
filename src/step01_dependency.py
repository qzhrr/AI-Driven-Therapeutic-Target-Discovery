"""
step01_dependency — Glioma-selective dependency (DepMap)

Cohen's d + Welch t + BH-FDR of glioma vs non-glioma CRISPR gene effect; emit 283 candidates.

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

# Glioma-selective dependency (DepMap)


import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Load model metadata
model = pd.read_csv(raw("Model.csv"))

# Load gene effect matrix
ge = pd.read_csv(raw("CRISPRGeneEffect.csv"), index_col=0)

# Clean gene names: "SYMBOL (Entrez)" -> SYMBOL
ge.columns = [c.split(" (")[0] for c in ge.columns]

# Get CRISPR model IDs
crispr_models = set(ge.index)

# Get glioma lines with CRISPR data
glioma = model[(model.OncotreeLineage=="CNS/Brain") & (model.OncotreePrimaryDisease=="Diffuse Glioma")].copy()
glioma["has_crispr"] = glioma.ModelID.isin(crispr_models)

glioma_ids = set(glioma[glioma.has_crispr].ModelID)
is_glioma = ge.index.isin(glioma_ids)

G = ge[is_glioma]
R = ge[~is_glioma]

n_g, n_r = len(G), len(R)
g_mean = G.mean()
r_mean = R.mean()
g_std = G.std()
r_std = R.std()

# Cohen's d
pooled_sd = np.sqrt(((n_g-1)*g_std**2 + (n_r-1)*r_std**2) / (n_g+n_r-2))
cohens_d = (g_mean - r_mean) / pooled_sd

# Welch's t-test
t, p = stats.ttest_ind(G, R, equal_var=False, nan_policy="omit")

res = pd.DataFrame({
    "gene": g_mean.index,
    "glioma_mean": g_mean.values,
    "rest_mean": r_mean.values,
    "diff": (g_mean - r_mean).values,
    "cohens_d": cohens_d.values,
    "t": t, "p": p,
}).set_index("gene")

# FDR correction
res = res.dropna(subset=["p"])
res["fdr"] = multipletests(res["p"], method="fdr_bh")[1]

# Common-essential exclusion list
ce = pd.read_csv(raw("CRISPRInferredCommonEssentials.csv"))
common_essentials = set(x.split(" (")[0] for x in ce.iloc[:,0])
res["is_common_essential"] = res.index.isin(common_essentials)

# Pan-essentiality feature
dep_frac_all = (ge < -0.5).mean(axis=0)
res["pan_ess_fraction"] = dep_frac_all.reindex(res.index).values

# Subset-dependency: fraction of GLIOMA lines where gene is a strong dependency
G_dep_frac = (G < -0.5).mean(axis=0)
res["glioma_dep_fraction"] = G_dep_frac.reindex(res.index).values

# Save checkpoint
res.reset_index().to_parquet(result("step01_dependency", "glioma_dependency_scored.parquet"))
print("checkpoint written:", res.shape)


import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Load data
model = pd.read_csv(raw("Model.csv"))

# Load gene effect matrix index to get CRISPR models
ge_index = pd.read_csv(raw("CRISPRGeneEffect.csv"), usecols=[0])
crispr_models = set(ge_index.iloc[:,0])

# Glioma lines with CRISPR data
glioma = model[(model.OncotreeLineage=="CNS/Brain") & (model.OncotreePrimaryDisease=="Diffuse Glioma")].copy()
glioma["has_crispr"] = glioma.ModelID.isin(crispr_models)

# Load full gene-effect matrix
ge = pd.read_csv(raw("CRISPRGeneEffect.csv"), index_col=0)
ge.columns = [c.split(" (")[0] for c in ge.columns]

# Cohort masks
glioma_ids = set(glioma[glioma.has_crispr].ModelID)
is_glioma = ge.index.isin(glioma_ids)

G = ge[is_glioma]
R = ge[~is_glioma]

# Per-gene stats
g_mean = G.mean()
r_mean = R.mean()
g_std = G.std()
r_std = R.std()
n_g, n_r = len(G), len(R)

# Pooled SD for Cohen's d
pooled_sd = np.sqrt(((n_g-1)*g_std**2 + (n_r-1)*r_std**2) / (n_g+n_r-2))
cohens_d = (g_mean - r_mean) / pooled_sd

# Welch's t-test
t, p = stats.ttest_ind(G, R, equal_var=False, nan_policy="omit")
res = pd.DataFrame({
    "gene": g_mean.index,
    "glioma_mean": g_mean.values,
    "rest_mean": r_mean.values,
    "diff": (g_mean - r_mean).values,
    "cohens_d": cohens_d.values,
    "t": t, "p": p,
}).set_index("gene")

# FDR correction
res = res.dropna(subset=["p"])
res["fdr"] = multipletests(res["p"], method="fdr_bh")[1]

# Common-essential exclusion list
ce = pd.read_csv(raw("CRISPRInferredCommonEssentials.csv"))
common_essentials = set(x.split(" (")[0] for x in ce.iloc[:,0])
res["is_common_essential"] = res.index.isin(common_essentials)

# Pan-essentiality as scored feature
dep_frac_all = (ge < -0.5).mean(axis=0)
res["pan_ess_fraction"] = dep_frac_all.reindex(res.index).values

# Subset-dependency: fraction of GLIOMA lines where gene is a strong dependency
G_dep_frac = (G < -0.5).mean(axis=0)
res["glioma_dep_fraction"] = G_dep_frac.reindex(res.index).values

# Revised candidate rule
def candidates_v2(d_floor=0.3, fdr_cut=0.05, mean_cut=-0.5, subset_cut=0.20):
    depish = (res.glioma_mean < mean_cut) | (res.glioma_dep_fraction >= subset_cut)
    m = (res.cohens_d < -d_floor) & (res.fdr < fdr_cut) & depish
    return res[m]

cand = candidates_v2(d_floor=0.3, fdr_cut=0.05, mean_cut=-0.5, subset_cut=0.20).copy()

# Add composite selectivity ranking score
cand["neglog_fdr"] = -np.log10(cand.fdr.clip(lower=1e-300))
cand["dependency_score"] = (-cand.cohens_d) * (-cand.glioma_mean.clip(upper=0)).clip(lower=0.01)
cand = cand.sort_values("cohens_d")

# Annotate archetype
cand["archetype"] = np.where(cand.glioma_mean < -0.5, "pan-glioma-dep",
                     np.where(cand.glioma_dep_fraction>=0.20, "subset-dep","selective-weak"))

# Save candidate table
cand_out = cand.copy()
cand_out.index.name = "gene"
cand_out.to_csv(result("step01_dependency", "glioma_candidates.csv"))


import matplotlib.pyplot as plt, numpy as np, pandas as pd, matplotlib as mpl

apply_figure_style()
res = pd.read_csv(result("step01_dependency", "depmap_glioma_dependency.csv"), index_col="gene")
cand = pd.read_csv(result("step01_dependency", "glioma_candidates.csv"), index_col="gene")
x = res["cohens_d"]; y = -np.log10(res["fdr"].clip(lower=1e-300))
fig, ax = plt.subplots(figsize=(7.2,5.6))
ax.scatter(x, y, s=5, c=META_GREY, alpha=0.30, linewidths=0, rasterized=True, label="all genes")
is_cand = res.index.isin(cand.index)
ax.scatter(x[is_cand], y[is_cand], s=14, c="#c0392b", alpha=0.75, linewidths=0,
           label=f"selective candidates (n={is_cand.sum()})")
CAND_C="#1f6f3d"; EXCL_C="#111111"
labels = {
 "JUN":(-1.45,8.6,CAND_C), "VRK1":(-1.5,6.6,CAND_C), "MET":(-1.32,2.6,CAND_C),
 "PTK2":(-0.98,5.7,CAND_C), "GPX4":(-1.2,4.1,CAND_C),
 "PDGFRA":(-1.0,0.9,CAND_C), "FGFR1":(-0.78,4.3,CAND_C),
 "CDC25B":(-0.33,5.3,CAND_C), "PTPN11":(-0.08,3.5,CAND_C),
 "SOX2":(-0.72,0.35,EXCL_C), "OLIG2":(-0.15,0.32,EXCL_C),
}
for g,(lx,ly,col) in labels.items():
    gx,gy = res.loc[g,"cohens_d"], -np.log10(max(res.loc[g,"fdr"],1e-300))
    incand = g in cand.index
    ax.scatter([gx],[gy], s=34, facecolors='none', edgecolors=(CAND_C if incand else EXCL_C),
               linewidths=1.4, zorder=6)
    ax.annotate(g, xy=(gx,gy), xytext=(lx,ly), fontsize=7, fontstyle="italic",
                fontweight=("bold" if not incand else "normal"), color=col, zorder=7,
                ha="center", arrowprops=dict(arrowstyle="-", lw=0.6, color=col, shrinkA=1, shrinkB=3))
from matplotlib.lines import Line2D
handles,_ = ax.get_legend_handles_labels()
handles += [Line2D([0],[0],marker='o',mfc='none',mec=CAND_C,ls='',ms=6,label='candidate (labeled)'),
            Line2D([0],[0],marker='o',mfc='none',mec=EXCL_C,ls='',ms=6,label='excluded control')]
ax.axvline(-0.3, ls="--", lw=0.8, color="#777")
ax.set_xlabel("Glioma selectivity  (Cohen's d;  \u2190 more glioma-selective)")
ax.set_ylabel("Significance  (\u2212log\u2081\u2080 FDR)")
ax.set_title("Glioma-selective CRISPR dependencies (DepMap 24Q4; 72 glioma vs 1,106 other lines)",
             fontsize=8, loc="left")
ax.legend(handles=handles, frameon=False, fontsize=6.5, loc="upper left")
ax.set_ylim(-0.6, None); ax.margins(0.04)
fig.tight_layout()
fig.savefig(result("step01_dependency", "dependency_volcano.png"), dpi=200)
