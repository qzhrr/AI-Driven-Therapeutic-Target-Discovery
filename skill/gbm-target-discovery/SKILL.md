---
name: gbm-target-discovery
description: >-
  Cancer small-molecule target-discovery pipeline that integrates genetic dependency
  (DepMap CRISPR), synthetic lethality, single-cell tumor-vs-microenvironment specificity,
  therapeutic-window safety, tumor enrichment, and druggability into a weighted, ranked
  candidate list — then uses an LLM ensemble to select category weights and nominate the
  strongest targets, with a gene-masked ablation that tests whether nominations are driven
  by evidence or by gene-name recognition. Use this skill whenever the user wants to
  discover, prioritize, rank, or nominate therapeutic targets (drug targets, dependencies,
  vulnerabilities, synthetic-lethal partners) for a cancer or other disease — especially
  when they mention DepMap, CRISPR dependency, tumor-vs-normal selectivity, tractability,
  a "therapeutic window", or want LLM-assisted target prioritization with bias controls.
  Also use it to reproduce or extend the glioblastoma (GBM) analysis in this repository,
  or to adapt the 12-step framework to a different cancer type or therapeutic lane.
license: MIT
---

# Cancer small-molecule target discovery

A reusable, evidence-integrating pipeline for nominating **novel, druggable targets** in a
cancer of interest. It was built and validated on glioblastoma (GBM), where it screened
17,916 genes down to 283 glioma-selective candidates, recovered 4 of 5 pre-declared positive
controls in the top 12 (PDGFRA #2, MET #4, VRK1 #7, FGFR1 #12), and surfaced novel,
evidence-carried leads (VRK1 and ELAVL1 survive gene-masking) plus a provocative
KIF2C–9p21/CDKN2A synthetic-lethality hypothesis.

This skill lets you (1) **reproduce or inspect** the GBM analysis, or (2) **adapt the
framework** to a new cancer type / therapeutic lane. The scientific value is the
*integration recipe* and the *bias-control ablations*, not any single dataset.

## When to use this skill

Trigger this skill when the task involves prioritizing or nominating therapeutic targets
from multi-omic evidence — e.g. "find drug targets in <cancer>", "rank CRISPR dependencies
by druggability and tumor selectivity", "which synthetic-lethal partners of <lesion> are
tractable", "reproduce the GBM target ranking", "adapt this pipeline to pancreatic cancer".
It is specifically useful when the user wants an *auditable, weighted* ranking with explicit
controls for LLM bias — not a single opaque score.

## The method (what makes it worth reusing)

Seven evidence categories are computed per gene, normalized to [0,1], and combined with
LLM-selected weights:

1. **dependency** — glioma/tumor-selective CRISPR gene effect (Cohen's *d*, Welch *t*, BH-FDR)
2. **selectivity** — penalty for pan-essential (housekeeping) genes
3. **synthetic_lethality** — genotype-stratified vulnerability to recurrent lesions + an
   expression-context paralog scan (catches partners silenced non-genetically, e.g. VRK1←VRK2)
4. **tumor_specificity** — single-cell malignant-vs-microenvironment expression
5. **tumor_enrichment** — tumor-vs-normal-tissue over-expression + prevalence
6. **therapeutic_window** — graded normal-tissue / normal-brain safety
7. **tractability** — small-molecule druggability (Open Targets / ChEMBL / DGIdb / UniProt),
   kept **separate** from a novelty axis so novelty never inflates the quality score

Two design choices are what make the output trustworthy and are the reusable core:

* **LLM weighting is blind to candidate scores.** Weights are assigned from category
  *definitions* only, firewalling the weighting from the ranking it produces. Run an
  expert vs non-expert framing to quantify how much domain framing (not data) drives weights.
* **Gene-masked nomination ablation.** Re-run the LLM nomination on the *same evidence* with
  gene names replaced by opaque codes. Candidates whose nomination frequency survives masking
  are evidence-carried; those that collapse were riding name recognition. This is the single
  most important validation and should be applied to any LLM-in-the-loop nomination.

`helpers` (auto-loaded from `kernel.py`) implements the reusable scoring primitives:
`cohens_d`, `benjamini_hochberg`, `minmax01`, `graded_percentile_safety`,
`weighted_composite`, and `masked_ablation_effect`. Full method: `docs/METHODS.md`.

## Reproduce the GBM analysis

The pipeline lives in `src/` at the repository root (one level up from this skill). It runs
in two modes (see `src/common.py`):

* **cached (default)** — reproduces the published results offline from the committed
  `data/cache/` and `results/`; no downloads, no API key.
* **live** (`GBM_LIVE=1`) — re-downloads raw inputs and re-queries the connectors / LLM
  (needs `anthropic` + `ANTHROPIC_API_KEY`); stochastic, so it won't match bit-for-bit.

```bash
python src/run_all.py --list          # show the 12 steps
python src/run_all.py                 # full pipeline, cached
python src/run_all.py --only 09 10 11a  # weights, nomination, gene-masked ablation
```

Each step is standalone (`python src/stepNN_*.py`) and writes to `results/stepNN_*/`.
The master output is `results/step10_nomination/final_candidate_ranking.csv` (283 genes,
all 7 category scores + composite + nomination frequency + gene-masking columns).

## Adapt to a new cancer type or lane

The 12-step structure is disease-agnostic; the parts to change are the cohort definition and
the domain framing:

1. **Cohort (step 1).** Change the DepMap lineage/disease filter (in `src/step01_dependency.py`)
   from `OncotreeLineage == "CNS/Brain"` & `OncotreePrimaryDisease == "Diffuse Glioma"` to your
   cancer's Oncotree terms. Keep the candidate rule (Cohen's *d* floor, FDR, dependency fraction).
2. **Recurrent lesions (step 3).** Replace the 9 GBM lesions (CDKN2A/MTAP del, PTEN/RB1/NF1 loss,
   EGFR/PDGFRA amp, TP53 mut) with your cancer's recurrent driver events.
3. **Normal-tissue references (steps 4–6).** Swap the single-cell atlas, normal-tissue panel, and
   tumor-vs-normal comparison for tissue-matched references (the GBM run uses brain-specific
   safety; a different cancer needs its own).
4. **Domain framing (steps 9–11a).** Update the expert system prompt and category definitions in
   `src/step09_weight_selection.py` / `src/step10_nomination.py` to your cancer and lane. The
   weighting, nomination, and gene-masked ablation logic carry over unchanged.

The tractability, feature-matrix assembly, weighting, nomination, and all ablation logic are
cancer-independent and need no changes.

## Outputs and interpretation

* `final_candidate_ranking.csv` — the master ranked table.
* `gene_masked_ablation.csv` — per-gene `freq_named`, `freq_masked`, and an `effect` class
  (`evidence-carried (robust)`, `NAME-INFLATED`, `NAME-SUPPRESSED`). **Read this before trusting
  any nomination.**
* `weights_aggregated.csv` / `weights_comparison.csv` — the weight vector and the expert-vs-
  non-expert framing effect.

**These are computational hypotheses, not validated targets.** Positive-control recovery
(known drivers ranking highly) is evidence the method works; a novel top pick is a lead to
test experimentally, not a conclusion. Report the gene-masking result alongside any nomination.
