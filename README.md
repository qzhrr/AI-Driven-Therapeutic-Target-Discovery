# AI-Driven Therapeutic Target Discovery

*A 12-step small-molecule target-discovery pipeline, developed and validated on glioblastoma (GBM).*

A 12-step computational pipeline that nominates **novel, druggable small-molecule
inhibition targets in glioblastoma (GBM)** by integrating a genetic-dependency core
(DepMap CRISPR) with a therapeutic-window filter (single-cell tumor-vs-microenvironment,
normal-tissue safety, tumor enrichment, tractability), then applying LLM-selected
category weights and an LLM nomination ensemble — with built-in ablations that test
whether the nominations are driven by *evidence* or by *gene-name recognition*.

The pipeline screens **17,916 genes → 283 glioma-selective candidates → a ranked
shortlist**. Its distinctive, deliberately provocative lead is **KIF2C**, whose
glioma dependency tracks the recurrent **9p21/CDKN2A** deletion as a synthetic-lethal
relationship — a novel, testable, and as-yet-unvalidated hypothesis.

> Single model throughout (no cross-model consensus). Feature definitions were
> specified by the scientist; the pipeline executes scoring, weighting, and nomination.
> **These are computational hypotheses, not validated targets** — see
> [Limitations](#limitations).

---

## Headline results

| | |
|---|---|
| Genes screened (DepMap 24Q4) | 17,916 |
| Glioma-selective candidates | **283** (236 pan-glioma-dependent, 47 subtype-dependent) |
| Top-ranked candidate | **KIF2C** (novel tier; 9p21/CDKN2A synthetic-lethality lead) |
| Robust novel picks (survive gene-masking) | **VRK1** (masked nom. freq 0.99), **ELAVL1** (0.92) |
| Name-inflated (collapse when masked) | KIF18A (1.00 → 0.00), CDK7 (1.00 → 0.03) |
| Positive controls recovered | 4/5 in the top 12 — PDGFRA #2, MET #4, VRK1 #7, FGFR1 #12 (PTPN11 #74, correctly down-weighted as common-essential) |

The **gene-masked ablation** (step 11a) is the methodological centerpiece: re-running
the LLM nomination on the *same evidence* with gene identities replaced by opaque
codes separates candidates whose support is real from those riding name recognition.

See [`docs/FINAL_REPORT.md`](docs/FINAL_REPORT.md) for the full scientific write-up and
[`docs/METHODS.md`](docs/METHODS.md) for complete methods, versions, and thresholds.

---

## The 12-step pipeline

| Step | Script | What it does |
|------|--------|--------------|
| 0  | `step00_fetch_data.py`            | Download all external inputs (DepMap, GBmap, GTEx, Lake, Xena) |
| 1  | `step01_dependency.py`            | Glioma-selective dependency: Cohen's *d* + Welch *t* + BH-FDR over CRISPR gene effect → 283 candidates |
| 3  | `step03_synthetic_lethality.py`   | Genotype-stratified SL vs 9 recurrent GBM lesions + expression-context paralog scan |
| 4  | `step04_tumor_vs_tme.py`          | Pseudobulk GBmap single-cell: malignant-vs-microenvironment specificity |
| 5  | `step05_therapeutic_window.py`    | Graded safety over GTEx + brain subregions + Lake normal-brain neurons |
| 6  | `step06_tumor_enrichment.py`      | TCGA-GBM vs GTEx-brain log2FC + prevalence (UCSC Xena Toil recompute) |
| 7  | `step07_tractability.py`          | Open Targets / ChEMBL / DGIdb / UniProt → graded tractability + a separate novelty axis |
| 8  | `step08_feature_matrix.py`        | Assemble the 283 × 12 feature matrix into 7 evidence categories; normalize |
| 9  | `step09_weight_selection.py`      | LLM distributes 100 points across the 7 categories (100 expert + 100 non-expert runs) |
| 10 | `step10_nomination.py`            | LLM nominates the strongest candidates (100 runs); assemble the master 283-gene ranking |
| 11 | `step11_ablations.py`             | Ablations: composite-only vs +LLM, expert vs non-expert, weight jitter, control recovery |
| 11a| `step11a_gene_masked_ablation.py` | **Blinded** re-nomination on anonymized candidates — evidence vs name recognition |
| 11b| `step11b_kif2c.py`                | KIF2C 9p21/CDKN2A SL deep-dive: genotype-selective dependency, mechanism, novelty ledger |

(Step 2 was merged into step 1 during the original analysis; the numbering is kept for
traceability with `docs/METHODS.md`.)

Run order and dependencies are encoded in [`src/run_all.py`](src/run_all.py).

---

## Quick start

```bash
# 1. Environment (conda or pip)
conda env create -f environment.yml && conda activate gbm-target-discovery
#   or:  python -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt

# 2. Reproduce the published results offline (cached mode — no downloads, no API keys)
python src/run_all.py

# Inspect / run individual steps
python src/run_all.py --list
python src/step09_weight_selection.py
python src/run_all.py --only 09 10 11a
```

Every script is runnable standalone (`python src/stepNN_*.py`) and writes its outputs
to `results/stepNN_*/`.

### Two run modes

The pipeline has exactly two stages that are not pure computation — connector queries
(`host.mcp`) and LLM calls (`host.llm`). Both are **cache-first**:

* **cached (default)** — reproduces the published results from the committed
  `data/cache/` (connector responses) and `results/` (including the exact 100-run LLM
  outputs). No network, no API key. Fully offline and deterministic.
* **live** (`GBM_LIVE=1`) — re-downloads raw inputs and re-queries the live APIs. The
  LLM steps need the `anthropic` package and `ANTHROPIC_API_KEY`. Because those steps
  are stochastic, a live re-run draws fresh samples and will not reproduce the committed
  values bit-for-bit.

```bash
export ANTHROPIC_API_KEY=sk-...      # live LLM steps only
GBM_LIVE=1 python src/run_all.py
```

> **Portability shim.** The project was originally executed inside Claude Science,
> where a `host` object provided `host.llm` / `host.mcp`. [`src/adapters.py`](src/adapters.py)
> re-implements that surface against the public Anthropic API so the code runs unchanged
> outside the platform; connector-dependent steps fall back to the committed cache.

---

## Repository layout

```
AI-Driven-Therapeutic-Target-Discovery/
├── src/                      # 16 runnable scripts + run_all.py
│   ├── common.py             #   paths, cache-first run-mode helpers
│   ├── adapters.py           #   portable host.llm / host.mcp shim
│   ├── figstyle.py           #   shared figure styling
│   ├── step00_fetch_data.py … step11b_kif2c.py
│   └── run_all.py            #   orchestrator (--from / --only / --list)
├── data/
│   ├── raw/                  # large external inputs (git-ignored; fetched by step00)
│   └── cache/                # small connector/LLM responses (committed → offline repro)
├── results/                  # per-step outputs: tables, figures, step reports
│   └── step01_dependency/ … step11b_kif2c/
├── docs/                     # METHODS.md, FINAL_REPORT.md, summary, figure captions
├── skill/gbm-target-discovery/   # Claude Science skill wrapping the pipeline
├── requirements.txt · environment.yml · LICENSE · .gitignore
```

### Data provenance

| Source | Version / access | Committed? |
|--------|------------------|------------|
| DepMap | 24Q4 Public (figshare 27993248): CRISPR/Chronos gene effect, 72 glioma vs 1,106 non-glioma lines | raw (fetched) |
| GBmap single-cell | CELLxGENE Census core atlas, 338,564 cells | raw (fetched, ~8.1 GB) |
| GTEx | v8, 54 tissues, bulk median TPM | cache (committed) |
| Lake normal brain | GEO GSE97930 snDrop-seq (frontal/visual cortex + cerebellar hem) | raw (fetched) |
| TCGA-GBM | UCSC Xena Toil recompute, 153 GBM vs 1,152 GTEx brain | cache/raw |
| Open Targets, ChEMBL, DGIdb, UniProt | tractability, mechanism, drug-gene, protein class | cache (committed) |

Raw inputs total ~8.6 GB and are **not committed** (see [`data/raw/README.md`](data/raw/README.md));
`step00_fetch_data.py` downloads them from source. The ~3.2 MB of connector/LLM caches
**are** committed so the API-dependent steps reproduce offline.

---

## Limitations

* **Hypotheses, not validated targets.** Every nomination is computational. The KIF2C–9p21
  synthetic-lethality lead in particular is novel and **experimentally unvalidated** — it is
  offered as a testable hypothesis, not a conclusion.
* **Single model, no consensus.** Weights and nominations come from one model family; the
  expert/non-expert and gene-masked ablations quantify the framing/name-recognition effects
  but do not remove them.
* **Small-molecule lane only.** Tumor-suppressor and CAR-T / surface-antigen targets are out
  of scope by design; a gene flagged as an SL-partner of a lost suppressor is an annotation,
  not a small-molecule nomination.
* **LLM steps are stochastic.** Exact reproduction requires the committed cache; a live re-run
  will differ at the margin.

## Citation & data terms

The pipeline **code** is MIT-licensed ([`LICENSE`](LICENSE)). The external datasets carry
their own licenses and terms of use — cite the original sources (DepMap, GBmap/CELLxGENE,
GTEx, Lake et al. 2018, TCGA/Xena, Open Targets, ChEMBL, DGIdb, UniProt) when using their
data. Versions and accessions are in [`docs/METHODS.md`](docs/METHODS.md).
