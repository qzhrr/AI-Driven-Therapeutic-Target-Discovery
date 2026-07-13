"""
step00_fetch_data — Data acquisition (all external downloads)

Downloads the large external inputs the pipeline needs into data/raw/ and derives
two small committed cache files (Xena sample lists, GTEx median expression).

Sources (versions per docs/METHODS.md):
  * DepMap 24Q4 Public (figshare article 27993248): CRISPRGeneEffect, Model,
    CRISPRInferredCommonEssentials, OmicsAbsoluteCNGene, somatic-mutation matrices.
  * GBmap core single-cell atlas (CELLxGENE Census, 338,564 cells) as .h5ad.
  * Lake et al. 2018 normal-brain snDrop-seq (GEO GSE97930), 3 cortical regions.
  * UCSC Xena Toil TcgaTargetGtex phenotype (to derive GBM-tumor / GTEx-brain sample lists).
  * GTEx v8 median expression (via the platform ``expression`` connector; cache-first).

Run modes (see common.py):
  * cached (default): raw downloads are skipped if the file already exists; the
    committed GTEx cache in data/cache/ is used as-is (no connector calls).
  * live (GBM_LIVE=1): re-fetch raw inputs and re-query the GTEx connector.

The raw inputs total ~8.6 GB (core_gbmap.h5ad alone is ~8.1 GB); they are
git-ignored. This script only downloads what is missing.
"""
import os
import sys
import ssl
import json
import gzip
import subprocess
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (  # noqa: E402
    DATA_RAW, LIVE, use_cache, raw, cache, result, require, ensure_dirs,
)
from adapters import host  # noqa: E402  portable host.mcp shim

ensure_dirs()

# TLS context tolerant of the occasional stale cert on public data mirrors.
_CTX = ssl.create_default_context()
_CTX.check_hostname = False
_CTX.verify_mode = ssl.CERT_NONE
_UA = {"User-Agent": "python-urllib"}


def _download(url: str, dest, timeout: int = 300) -> None:
    """Fetch url -> dest unless the file already exists (skip in cached mode)."""
    dest = str(dest)
    if os.path.exists(dest) and not LIVE:
        print(f"  cached: {os.path.basename(dest)} ({os.path.getsize(dest):,} bytes)")
        return
    req = urllib.request.Request(url, headers=_UA)
    with urllib.request.urlopen(req, timeout=timeout, context=_CTX) as r:
        data = r.read()
    with open(dest, "wb") as f:
        f.write(data)
    print(f"  fetched: {os.path.basename(dest)} ({len(data):,} bytes)")


# ---------------------------------------------------------------------------
# 1. DepMap 24Q4 (figshare) — CRISPR gene effect, model metadata, essentials, CN, mutations
# ---------------------------------------------------------------------------
DEPMAP_FILES = {
    51064667: "CRISPRGeneEffect.csv",
    51065297: "Model.csv",
    51064916: "CRISPRInferredCommonEssentials.csv",
    51065303: "OmicsAbsoluteCNGene.csv",
    51065747: "OmicsSomaticMutationsMatrixDamaging.csv",
    51065750: "OmicsSomaticMutationsMatrixHotspot.csv",
}
print("[1/5] DepMap 24Q4 (figshare)")
for fid, name in DEPMAP_FILES.items():
    _download(f"https://ndownloader.figshare.com/files/{fid}", raw(name), timeout=600)

# ---------------------------------------------------------------------------
# 2. GBmap core single-cell atlas (CELLxGENE) — ~8.1 GB
# ---------------------------------------------------------------------------
print("[2/5] GBmap core atlas (CELLxGENE)")
_gbmap = raw("core_gbmap.h5ad")
if os.path.exists(_gbmap) and not LIVE:
    print(f"  cached: core_gbmap.h5ad ({os.path.getsize(_gbmap):,} bytes)")
else:
    # curl streams the large h5ad more reliably than urllib.
    subprocess.run([
        "curl", "-sL", "-m", "3600", "-o", str(_gbmap),
        "-w", "  HTTP %{http_code} size %{size_download} time %{time_total}s\n",
        "https://datasets.cellxgene.cziscience.com/861acfd8-25f0-418b-a445-aa96da232827.h5ad",
    ], check=True)

# ---------------------------------------------------------------------------
# 3. Lake et al. 2018 normal-brain snDrop-seq (GEO GSE97930) — 3 regions
# ---------------------------------------------------------------------------
print("[3/5] Lake 2018 normal brain (GSE97930)")
_GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97930/suppl"
for region in ("FrontalCortex", "VisualCortex", "CerebellarHem"):
    fn = f"GSE97930_{region}_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz"
    _download(f"{_GEO}/{fn}", raw(fn), timeout=600)

# ---------------------------------------------------------------------------
# 4. UCSC Xena Toil phenotype -> GBM-tumor / GTEx-brain sample lists (committed cache)
# ---------------------------------------------------------------------------
print("[4/5] Xena TcgaTargetGtex phenotype -> sample lists")
if use_cache(cache("xena_samples.json")):
    print("  cached: xena_samples.json")
else:
    import pandas as pd
    pheno_gz = raw("TcgaTargetGTEX_phenotype.txt.gz")
    _download("https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz",
              pheno_gz, timeout=300)
    ph = pd.read_csv(pheno_gz, sep="\t", encoding="latin-1")
    disease, stype, study = "primary disease or tissue", "_sample_type", "_study"
    gbm_tumor = ph[(ph[disease] == "Glioblastoma Multiforme") &
                   (ph[stype] == "Primary Tumor")]["sample"].tolist()
    gtex_brain = ph[(ph[study] == "GTEX") &
                    (ph[disease].str.startswith("Brain", na=False))]["sample"].tolist()
    json.dump({"gbm_tumor": gbm_tumor, "gtex_brain": gtex_brain},
              open(cache("xena_samples.json"), "w"))
    print(f"  GBM tumor: {len(gbm_tumor)}  GTEx brain: {len(gtex_brain)}")

# ---------------------------------------------------------------------------
# 5. GTEx v8 median expression via the platform 'expression' connector (cache-first)
# ---------------------------------------------------------------------------
# The connector query is candidate-scoped, so it depends on step01's candidate
# list. In cached mode the committed gtex_median*.json are used as-is; this block
# only fires in live mode (GBM_LIVE=1) and requires the platform connector.
print("[5/5] GTEx v8 median expression (connector, cache-first)")
if use_cache(cache("gtex_median.json")) and use_cache(cache("gtex_median_alias.json")):
    print("  cached: gtex_median.json, gtex_median_alias.json")
else:
    import csv
    require(result("step01_dependency", "glioma_candidates.csv"),
            "Run step01_dependency.py first (GTEx query is scoped to its candidates).")
    require(cache("gtex_resolved.json"),
            "gtex_resolved.json (committed) maps candidate symbols to GTEx gencode IDs.")

    with open(result("step01_dependency", "glioma_candidates.csv")) as f:
        cand_syms = [r["gene"] for r in csv.DictReader(f)]

    res = json.load(open(cache("gtex_resolved.json")))
    genes = res["genes"]
    gencode_ids = [g["gencodeId"] for g in genes]
    sym_by_gencode = {g["gencodeId"]: g["geneSymbol"] for g in genes}
    resolved_syms = {g["geneSymbol"] for g in genes}
    unresolved = [s for s in cand_syms if s not in resolved_syms]
    print("  unresolved (no GTEx v8 id):", unresolved)

    me = host.mcp("expression", "gtex_median_expression", gencode_ids=gencode_ids)
    rows = me.get("medianGeneExpression", me.get("data", me)) if isinstance(me, dict) else me
    json.dump({"rows": rows, "sym_by_gencode": sym_by_gencode, "unresolved": unresolved},
              open(cache("gtex_median.json"), "w"), default=str)

    # Recover legacy-symbol aliases GTEx v8 knows under an older name.
    alias = {"POLR1H": "ZNRD1", "ODR4": "C1orf27", "EPRS1": "EPRS", "POLR1G": "CD3EAP",
             "HJV": "HFE2", "NARS1": "NARS", "H3C8": "HIST1H3H", "TARS1": "TARS",
             "NOPCHAP1": "C12orf45"}
    res2 = host.mcp("expression", "gtex_resolve_genes", gene_ids=list(alias.values()))
    recs2 = res2.get("genes", [])
    rev = {v: k for k, v in alias.items()}
    gencode2 = [r["gencodeId"] for r in recs2]
    sym2 = {r["gencodeId"]: rev.get(r["geneSymbol"], r["geneSymbol"]) for r in recs2}
    if gencode2:
        me2 = host.mcp("expression", "gtex_median_expression", gencode_ids=gencode2)
        rows2 = me2["medians"]
        json.dump({"rows": rows2, "sym_by_gencode": sym2},
                  open(cache("gtex_median_alias.json"), "w"), default=str)

print("\nData acquisition complete.")
