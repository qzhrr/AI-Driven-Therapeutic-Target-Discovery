# GBM Target Discovery — Summary (v4)

*Generated from FINAL_REPORT.md (v4), METHODS.md, final_candidate_ranking.csv (283 candidates), and
gene_masked_ablation.csv. All statistics recomputed from raw source tables; independently re-derived in R_final_review.md.*

---

## 1. Hook

Starting from a whole-genome CRISPR screen, an LLM-weighted, dependency-first pipeline does not just rank glioblastoma
targets — it **triages** them by *kind of evidence*, using a blinded gene-masking test that removes the model's own
name-recognition bias and, in doing so, demotes its own top-ranked novel pick.

## 2. Structured abstract  (<200 words)

**Background.** Glioblastoma (GBM) has a ~15-month median survival and few novel druggable targets; naive scoring
pipelines tend to rediscover already-drugged oncogenes and reward highly-expressed genes rather than genetic
dependencies.

**Approach.** We ranked 17,916 genes from a glioma CRISPR dependency screen (DepMap 24Q4) down to 283
glioma-selective candidates, scored them on seven evidence categories (weighted by a domain-expert LLM) plus an
expression-context synthetic-lethality (SL) layer, and used a large language model to nominate targets across 100
runs. A blinded **gene-masking ablation** re-ran nomination on anonymised candidates to separate evidence from name
recognition.

**Results.** Four of five declared dependency positive controls recovered in the top 12 (Mann–Whitney p = 3.1×10⁻⁵).
Masking removed the model's name bias: KIF18A collapsed from a named nomination frequency of 1.00 to 0.00, while VRK1
held at 0.99. Three leads survive as a portfolio — **KIF2C** (a novel, unvalidated 9p21/CDKN2A-deletion synthetic
lethality), **VRK1** (a de-risked, still-undrugged small-molecule pick), and **ELAVL1/HuR** (evidence-carried, a
targeted-degrader opportunity).

**Conclusion.** The framework's contribution is a *calibrated triage*, not a target list — and it is self-skeptical
enough to remove its own name-inflated candidate.

## 3. Executive summary

**The problem.** GBM remains near-uniformly lethal and its standard of care has changed little in two decades. Two
failure modes plague in-silico target scoring: (i) it rewards genes that are merely highly expressed rather than
genetically required, and (ii) it rediscovers genes that are famous because they are already drugged. This pipeline
is built to avoid both.

**The pipeline in one breath.** A whole-genome CRISPR dependency screen (17,916 genes, DepMap 24Q4) is filtered to
283 glioma-selective dependencies, each scored across **seven categories** — dependency, selectivity, synthetic
lethality, tumour specificity, tumour enrichment, therapeutic window, tractability — combined into a composite whose
category weights are set *blind* by a domain-expert LLM. An expression-context SL layer (SL-v2) adds paralog/
co-complex partners, and an LLM nominates targets over 100 runs. Category weights sum to 1 and the composite
reconstructs from weights × categories to floating-point identity (max abs diff 1.1×10⁻¹⁶) — the ranking is exactly
the stated function of the stated inputs, with no post-hoc fudge factor.

**Headline: a three-way triage, not a leaderboard.**

| Lead | rank | named→masked nom. | tractability | kind of lead | biology |
|------|:----:|:-----------------:|:------------:|--------------|---------|
| **KIF2C** (MCAK) | 1 | 0.86 → 0.29 | 0.60 | **novel, unvalidated** | mitotic kinesin; selectively required in 9p21/CDKN2A-deleted glioma (d = −0.70, MWU p = 8.3×10⁻³, n = 44/28) — a synthetic-lethal pairing with no prior art |
| **VRK1** | 7 | 1.00 → 0.99 | 0.80 | **de-risked** | serine/threonine kinase; a dependency positive control *and* undrugged — the cleanest novel small-molecule pick; SL with paralog VRK2 |
| **ELAVL1** (HuR) | 22 | 0.98 → 0.92 | 0.55 | **chemically occupied** | RNA-binding protein; strong evidence but low classical tractability — a targeted-degrader opportunity, not a classical inhibitor |

VRK1 is simultaneously a positive control and a novelty-tier (undrugged) gene: it is the *bridge* case — known
biology that de-risks a still-open drug-discovery opportunity. It is deliberately **not** presented as a brand-new
discovery.

**The blinded test (the distinctive result).** Because scores were byte-identical between the named and masked arms,
only the gene's name changed. Masked nomination frequency tracks the pipeline's own selectivity evidence far better
than named frequency (Spearman ρ = 0.12 named → 0.57 masked). The test demoted the pipeline's own headline: KIF18A,
nominated in 100/100 named runs, fell to 0/100 when masked — it was nominated on name recognition, not evidence, and
was removed from the portfolio. KIF2C itself drops (0.86 → 0.29), which is why it is labelled *novel-but-unvalidated*
rather than a confirmed hit.

**What the LLM layer added — and its measured bias.** Beyond re-ranking, masking exposed the model's priors as
*data*: CDK7 (named 1.00 → masked 0.03) was a name-recognition artefact, not a hidden gem; JUN was name-vetoed
(0.03) then nominated unanimously (1.00) once anonymised; E2F3 was flagged "degrade, don't inhibit"; PDCD5 was
flagged as a tumour-suppressor reactivation target. The LLM adds directionality and modality judgement the
deterministic composite cannot, but it also imports name bias — which the masking control makes visible and
correctable.

**Expert vs non-expert.** Two independently-elicited weight vectors rank candidates similarly overall (Spearman
ρ = 0.917) but diverge by design: the expert weighting favours druggability and selectivity; the non-expert favours
raw proliferation/expression. The split is the quantitative statement of "domain expertise shifts the ranking toward
tractable, selective targets."

**Robustness.** Weight jitter of ±0.03 leaves the ranking essentially unchanged (Spearman ρ = 0.975; 89% of the
top-20 retained). Leave-one-out identifies selectivity and tractability as the load-bearing categories. The SL-v2
upgrade preserves the whole ranking (ρ = 0.993) while lifting VRK1 from #9 to #7. All 1,497 LLM nominations fall
within the 283-candidate universe — zero hallucinated genes.

**Caveats (faithful to the sources, incl. independent QC).**
- **KIF2C is a computationally-derived, unvalidated hypothesis.** The 9p21-SL pairing needs isogenic CDKN2A
  knockout/rescue confirmation.
- **The 9p21 attribution is block-level, not CDKN2A-specific.** CDKN2A and MTAP deletion are collinear (r = 0.97;
  one discordant line), so the biomarker is best stated as "9p21-block loss," and MTAP-synthetic-lethality (PRMT5/
  MAT2A) is the mechanistically nearest prior art.
- **Masking is a lower bound on name influence.** The score profile (tractability, novelty) can still leak identity;
  for the novel picks the leak is small (VRK1 tractability 0.80, ELAVL1 0.55 — not maxed), but the "evidence-only"
  reading is not pure.
- **Nomination is nested in the composite's own top-100**, so it re-ranks within the composite rather than
  supplying a fully independent signal (weight selection was firewalled, so this is not circular).
- **Masked-arm dropout** (87/100 valid runs vs 100/100 named) is itself non-random — the model refused more often
  when it could not anchor on names.
- 2D-culture CRISPR dependencies; IDH-wt-skewed cohort (71/72 IDH-wt, 51 GBM of 72 lines); BBB penetrance deferred
  to chemistry; single model (Claude), no cross-model consensus.
- Independent QC verdict: **GO-WITH-FIXES** — arithmetic reproduces from raw data; the one required fix (the swapped
  positive-control set) is applied here (4/5, p = 3.1×10⁻⁵, replacing the retracted 5/5, p = 3.5×10⁻⁸).

## 4. Methods at a glance

| Step | Data source (version) | Output |
|------|-----------------------|--------|
| Dependency screen | DepMap 24Q4 CRISPR (Chronos) | glioma-selective Cohen's d, FDR (17,916 genes) |
| Candidate filter | d ≤ −0.30 & FDR < 0.05 | 283 candidates |
| Expression / enrichment | GBmap, GTEx v8, TCGA-GBM | tumour specificity & enrichment scores |
| Synthetic lethality | genotype SL + SL-v2 expression-context paralog/co-complex | SL category; VRK1←VRK2, #9→#7 |
| Tractability / novelty | Open Targets, ChEMBL, DGIdb, UniProt | tractability & novelty-tier |
| Composite weighting | domain-expert LLM (blind to scores), 7 categories, Σweights=1 | expert_composite, rank_all |
| Nomination | LLM, 100 runs, top-100 by composite | nom_frequency (0 hallucinations / 1,497) |
| **Step 11a — gene-masked ablation** | anonymised candidates, 100 runs | freq_named vs freq_masked (87/100 valid) |

*Source of truth: raw tables (final_candidate_ranking.csv, gene_masked_ablation.csv) > FINAL_REPORT.md (v4) /
METHODS.md > derived files. Numbers that appear in the older demo script and README (5/5 controls; p = 4.3×10⁻⁸ or
3.5×10⁻⁸; KIF18A as a novel target) are superseded and are not used here.*
