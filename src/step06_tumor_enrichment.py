"""
step06_tumor_enrichment — Tumor expression enrichment (TCGA-GBM)

TCGA-GBM vs GTEx-brain log2FC + prevalence via UCSC Xena Toil recompute.

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

# Tumor expression enrichment (TCGA-GBM via Xena)


import pandas as pd, numpy as np, json
import urllib.request, ssl
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# Load candidates
cand = pd.read_csv(result("step01_dependency", "glioma_candidates.csv"))["gene"].tolist()

# Download probemap
ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
req = urllib.request.Request("https://toil.xenahubs.net/download/probeMap%2Fgencode.v23.annotation.gene.probemap",
                              headers={"User-Agent": "python-urllib"})
with urllib.request.urlopen(req, timeout=120, context=ctx) as r:
    open("/tmp/gtex_probemap", "wb").write(r.read())

pm = pd.read_csv("/tmp/gtex_probemap", sep="\t")

alias = {"EPRS1": "EPRS", "H3C8": "HIST1H3H", "HJV": "HFE2", "NARS1": "NARS", "NOPCHAP1": "C12orf45",
         "ODR4": "C1orf27", "POLR1G": "CD3EAP", "POLR1H": "ZNRD1", "SGO1": "SGOL1", "SPOUT1": "C9orf114", "TARS1": "TARS"}
rev = {v: k for k, v in alias.items()}
extra = pm[pm["gene"].isin(alias.values())].copy()
extra["sym"] = extra["gene"].map(lambda g: rev.get(g, g))
cand_pm = pm[pm["gene"].isin(cand)].copy(); cand_pm["sym"] = cand_pm["gene"]
allmap = pd.concat([cand_pm[["id", "sym"]], extra[["id", "sym"]]]).drop_duplicates("id")
allmap.to_csv("/tmp/ensg_sym.csv", index=False)
with open("/tmp/target_ensg.txt", "w") as f:
    for i in allmap["id"]: f.write(i + "\n")

# Stream and filter the TPM matrix
import subprocess, gzip, io

target_ids = set(allmap["id"])

req2 = urllib.request.Request("https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz",
                               headers={"User-Agent": "python-urllib"})
with urllib.request.urlopen(req2, timeout=1800, context=ctx) as r:
    with gzip.open(r, "rt") as gz:
        lines = []
        header = gz.readline()
        lines.append(header)
        for line in gz:
            gene_id = line.split("\t", 1)[0]
            if gene_id in target_ids:
                lines.append(line)

with open("/tmp/gbm_candidate_tpm.tsv", "w") as f:
    f.writelines(lines)

# Load filtered matrix
mat = pd.read_csv("/tmp/gbm_candidate_tpm.tsv", sep="\t", index_col=0)
ensg_sym = pd.read_csv("/tmp/ensg_sym.csv").set_index("id")["sym"].to_dict()
mat.index = [ensg_sym.get(i, i) for i in mat.index]
mat = mat[~mat.index.duplicated(keep="first")]

samp = json.load(open(cache("xena_samples.json")))
gbm_t = [s for s in samp["gbm_tumor"] if s in mat.columns]
gtx_b = [s for s in samp["gtex_brain"] if s in mat.columns]

T = mat[gbm_t]; N = mat[gtx_b]

feat = pd.DataFrame(index=mat.index)
feat["tcga_gbm_median_log2tpm"] = T.median(axis=1)
feat["gtex_brain_median_log2tpm"] = N.median(axis=1)
feat["tcga_vs_normal_log2fc"] = T.median(axis=1) - N.median(axis=1)
n75 = N.quantile(0.75, axis=1)
feat["tcga_pct_over_normal"] = (T.gt(n75, axis=0)).mean(axis=1)

pvals = []
for g in mat.index:
    try:
        pvals.append(mannwhitneyu(T.loc[g].dropna(), N.loc[g].dropna(), alternative="two-sided").pvalue)
    except Exception:
        pvals.append(np.nan)
feat["tcga_vs_normal_p"] = pvals
ok = feat["tcga_vs_normal_p"].notna()
feat.loc[ok, "tcga_vs_normal_fdr"] = multipletests(feat.loc[ok, "tcga_vs_normal_p"], method="fdr_bh")[1]

feat.to_csv(result("step06_tumor_enrichment", "tcga_tumor_enrichment.csv"))
