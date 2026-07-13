# METHODS — GBM small-molecule target discovery pipeline
*Cross-dataset integration + single-model LLM weighting/nomination. All computation by Claude for Science (single model; no cross-model consensus).*

## 0. Design
Goal: nominate **novel, druggable small-molecule inhibition targets** in glioblastoma (GBM) by integrating a
genetic-dependency core (DepMap) with a therapeutic-window filter (normal-tissue, tumor, single-cell, tractability),
then applying LLM-selected category weights and LLM nomination. Lane = small molecule (no CAR-T). Cancer = GBM.
Single model throughout. The scientist owned all feature definitions; the agent executed scoring, weighting, nomination.

## 1. Data sources & versions
| Source | Version / access | Use |
|---|---|---|
| DepMap | 24Q4 Public (figshare article 27993248, 2024-12-10); CRISPR/Chronos gene effect; 72 glioma vs 1,106 non-glioma lines |
| GBmap single-cell | CELLxGENE Census core atlas 338,564 cells (collection 999f2a15-3d7e-440b-96ae-2c806799c08c; DOI 10.1101/2022.08.27.505439) |
| GTEx | v8, 54 tissues, bulk median TPM (GTEx portal API v2) |
| Lake normal brain | GSE97930 snDrop-seq (Lake et al. 2018; frontal/visual cortex + cerebellar hem) |
| TCGA-GBM tumor | UCSC Xena Toil TcgaTargetGtex RSEM gene TPM (153 GBM tumor vs 1,152 GTEx brain, uniform recompute) |
| Open Targets | GraphQL (small-molecule tractability buckets + known-drug counts) |
| ChEMBL | EBI API (mechanism action_types + IC50/Ki/Kd/EC50 pChEMBL>=6) |
| DGIdb | v5 GraphQL (drug-gene interactions) |
| UniProt | REST (protein class + subcellular localization) |

## 2. Step 1 — Glioma-selective dependency (DepMap)
- Cohort: **all glioma** for power, subtype visible in metadata — 72 glioma CRISPR lines (51 GBM) vs 1,106 non-glioma;
  matrix 1,178 × 17,916. CRISPR/Chronos only (RNAi skipped per decision).
- Per gene: Cohen's d (glioma vs rest) + Welch t + BH-FDR. Pan-essentiality **scored** (`pan_ess_fraction`), NOT
  hard-filtered.
- Candidate rule: d_floor 0.3, subset_cut 0.20, FDR<0.05, dependency = (glioma_mean < −0.5 OR glioma_dep_fraction ≥ 0.20).
  → **283 candidates** (236 pan-glioma-dep, 47 subset-dep; 179 common-essential-flagged, 104 not).
- Controls: VRK1/MET/PDGFRA/FGFR1 retained; EGFR correctly excluded (CRISPR de-amplification artifact, rank 17,834/17,916);
  TSGs/master-TFs out.

## 3. Step 3 — Synthetic-lethality / genotype-stratified layer
- 9 recurrent GBM lesions (CDKN2A/CDKN2B/MTAP deletion, PTEN/RB1/NF1 loss, EGFR/PDGFRA amplification, TP53 mutation);
  CN thresholds del<0.5, amp>=6 copies.
- Pan-cancer association detection + glioma-direction-consistency + glioma-prevalence anchor (floor **8%** per decision).
- 67 candidates carry an anchored SL signal; 48 are SL-partners of a lost tumor-suppressor (annotation/bonus only —
  tumor suppressors themselves are OFF-lane for small-molecule inhibition).
- **v2 — expression-context SL detector (paralog/co-complex).** The genotype module above tests only copy-number/
  mutation lesions, so it is blind to partners silenced epigenetically/transcriptionally (e.g. VRK1←VRK2). Added an
  expression-context scan: correlate each candidate's CRISPR profile against genome-wide expression; a paralog whose
  low expression predicts stronger dependency is a candidate SL partner (44 significant hits). The 10 symbol-root paralog
  hits were verified against **Ensembl Compara** (`ensembl_homology`, paralogues); 9 are significant true paralogs
  (ACTL6A←ACTL6B drops on significance). **CHMP3←CHMP2A is retained** as an obligate ESCRT-III co-complex co-dependency
  (not buffering) per user decision. Integration: `sl_effect = max(genotype, expression-paralog d)`, **winsorized at the
  strongest genotype SL (E2F3, d=1.53)** before min-max normalization. Whole-ranking Spearman vs genotype-only = **0.993**;
  6 of the 9 credible paralog SLs were previously scored 0.0. VRK1 rose #9→#7; CHMP3 #255→#213.
  *Column note:* the winsorization is applied to the composite-feeding feature `f_sl_strength` (all SLs with |d|≥1.53 saturate
  to 1.0); the raw display column `sl_effect` retains the pre-cap magnitudes (FAM50A 3.14, CDS2 2.31, …) and does not enter
  the composite. This is a display/feature distinction with no scoring impact.

## 4. Step 4 — Tumor-vs-TME deconvolution (GBmap single-cell)
- GBmap core atlas (338,564 cells). Pseudobulk per cell-type compartment; malignant vs microenvironment.
- Per candidate: `sc_malignant_specificity`, `sc_tme_dominated_flag`. Compartment label renamed to
  **tumor-cell-expressed | TME-dominated** (expression compartment, not function). 43 TME-dominated, 240 tumor-expressed,
  5 not in atlas (imputed neutral downstream).
- Treatment (per decision): TME-dominated penalty + malignant-specificity as a modest bonus.

## 5. Step 5 — Therapeutic window (safety)
- GTEx v8 body-wide + 13 brain subregions; GBmap normal compartments; Lake GSE97930 normal-brain neurons.
- Graded percentile safety (rank-normalized [0,1], higher = safer), NOT an absolute-TPM cutoff.
  `safety_composite = mean(safety_vital_organ, safety_brain_bulk, safety_brain_neuron)`; body-wide max reported, excluded
  from the composite (equal-weight over 3 axes, locked by decision).

## 6. Step 6 — Tumor expression enrichment (TCGA-GBM)
- UCSC Xena Toil TcgaTargetGtex (uniform RSEM recompute → no batch confound): 153 GBM tumor vs 1,152 GTEx brain.
- `tcga_vs_normal_log2fc`, `tcga_pct_over_normal`, Welch p + BH-FDR.

## 7. Step 7 — Tractability + druggability (two independent axes)
- Open Targets SM tractability buckets; ChEMBL mechanism + genuine potent-compound count; DGIdb interactions; UniProt class/loc.
- **`tractability_score`** (0–1 graded ladder): HQ pocket/ligand +0.35, structure-with-ligand +0.25, druggable family +0.20,
  ≥10 genuine potent compounds +0.20.
- **`novelty_score`** kept a SEPARATE axis (approved drug does NOT down-weight; it only lowers novelty).
- **ChEMBL data-quality fix:** raw compound count inflated screening-only targets (LMNA 6,192 → 1 on genuine
  IC50/Ki/Kd/EC50; real drug targets unchanged: MET 6,729, FGFR1 5,414, CDK2 2,409).

## 8. Step 8 — Feature matrix (283 × 12 scored features → 7 categories)
- Categories: dependency, selectivity (pan-ess penalty), synthetic-lethality (bonus), tumor-specificity, tumor-enrichment,
  therapeutic-window, tractability. **Novelty is a separate reported axis, not one of the 7 weighted categories.**
- Normalization: rank-percentile (unbounded), keep-as-is ([0,1] features), min-max zero-preserving (SL). Penalties inverted.
- 5 GBmap-missing genes imputed neutral (0.5). Category score = mean of member features.
- Correlation check confirmed category structure absorbs redundancy (sl_strength~sl_breadth ρ=0.98) and preserves
  designed tensions (dep_breadth~pan_ess ρ=−0.96; tractability~novelty ρ=−0.53).

## 9. Step 9 — LLM weight selection (single-model ensemble)
- Model: **Claude Opus (`claude-opus-4-8`, session model)** — literal "just you," no consensus.
- Task: distribute 100 points across the 7 categories; forced `tool_use` schema; **blind to candidate scores**
  (weights on domain reasoning alone).
- Two conditions × 100 runs: **expert** (GBM SM target-discovery framing + domain category definitions) vs
  **non-expert** (generic feature framing, domain cues stripped).
- Constraint: Opus deprecates the `temperature` parameter → default sampling (documented; across-category mean CV 0.055).
  Per-category, the distribution is non-degenerate for six of seven categories, but **dependency returned exactly 24 points in
  all 100 expert runs (std=0, CV=0)** — so the "100-run ensemble" is effectively a single value for the highest-nominal-weight
  category. Low practical impact (dependency drives only ~6% of ranking influence), but the degeneracy is disclosed rather than
  masked by the aggregate CV.
- **Final weight vector (expert median, Σ=1):** dependency 0.24, selectivity 0.18, tractability 0.16, synthetic-lethality 0.12,
  therapeutic-window 0.12, tumor-specificity 0.11, tumor-enrichment 0.07.

## 10. Step 10 — LLM target nomination
- **Novelty correction adopted (from B's review, independently verified):** the Step-7/8 novelty label mislabeled 57 genes as
  approved-drugged via `ot_sm_approved_drug OR dgidb_n_approved_drugs>0`; only **21** carry the authoritative Open Targets
  approved-SM flag (the other 36 were DGIdb "approved-drug" records on undrugged genes — JUN, PTPN11, TEAD1, VCP, RAC1, GRB2,
  PTK2 — where a promiscuous approved drug creates a spurious link). Replaced with a 3-tier axis: **21 approved_sm_drug /
  15 has_sm_chemistry / 247 novel.** This changes no composite score (novelty was always separate) — only tier labels.
- **Posture P2-soft (user decision):** all 283 ranked in one list, tiers as visible tags, NO gating — faithful to the original
  "report, don't down-weight" ruling.
- Nomination: Opus reviewed the top 100 by expert composite each run, nominated ~15 with structured reason + **modality**
  (inhibit/degrade/reactivate/unknown) + confidence. **100 runs → nomination frequency.** 1,497 nominations, 0 hallucinations, 0 errors.

## 11. Step 11 — Ablations + validation
composite-only vs +LLM (ρ=0.49/0.55); expert vs non-expert weights (ρ=0.917, druggable vs proliferation split);
weight jitter ±0.03 (ρ=0.975, 89% top-20 retention); leave-one-out (selectivity + tractability load-bearing);
positive-control recovery on the pre-declared control set {PDGFRA, MET, VRK1, FGFR1, PTPN11}: 4/5 in top-12 (PTPN11 at #74,
correctly down-weighted as common-essential), MWU p=3.1e-5; EGFR/CDK4/CDK6 correctly excluded upstream).

### 11a. Gene-masked nomination ablation (blinded test)
- **Design:** two arms, identical harness (`claude-opus-4-8`, same system/tool/prompt, 100 runs each), differing in ONE
  variable — in the **masked** arm each candidate's gene symbol is replaced by an opaque code (e.g. `CX657`) while every
  score, category value, and pipeline flag is byte-identical. Named arm reproduced the original Step-10 consensus core
  exactly (harness is faithful). Valid runs: named 100/100, masked 87/100 (frequencies normalized per valid runs).
- **Result:** masked frequency tracks the pipeline's own evidence *better* than named (selectivity ρ 0.12→0.57; composite
  0.47→0.72); named-vs-masked correlate only ρ=0.32. Three behavior classes: **evidence-carried** (Δ≈0: VRK1, MET, PDGFRA,
  CDK2, ELAVL1, EGLN1, E2F3, KIF18B, FGFR1), **name-inflated** (drop when masked: KIF18A −1.00, CDK7 −0.97, PTPN11 −0.60,
  KIF2C −0.57, TEAD1 −0.53, WDR77 −0.40), **name-suppressed** (rise when masked: JUN +0.97, AMY2A +0.82, FERMT2 +0.53).
- **Effect on headline:** VRK1 and ELAVL1 are the evidence-carried novel nominations; KIF2C is borderline; **KIF18A was
  carried by name recognition and is removed from the novel headline.** Columns `freq_masked`, `masked_delta`,
  `masked_effect_class` added to the master ranking.

### 11b. KIF2C — 9p21/CDKN2A-deletion synthetic-lethality + novelty assessment
- **Genotype-selective dependency.** KIF2C dependency is enriched in 9p21.3-deleted glioma lines: CDKN2A-del mean gene-effect
  −0.503 vs −0.353 intact (**glioma Cohen's d=−0.70, Mann-Whitney p=8.3e-3**, n=44 del / 28 intact; 45% vs 14% strongly
  dependent; Welch t on the same split p=2.6e-3). Replicates as a pan-cancer anchored SL (**d=−0.46, raw p=2.8e-10,
  BH-FDR=7.9e-8, n=298**) and extends to co-deleted CDKN2B (d=−0.41) and MTAP
  (d=−0.37) — the full 9p21 block. Cohort: 61% CDKN2A-deleted, all 44 IDH-wt.
- **Mechanism.** Not proliferation-driven (glioma r=+0.21 n.s.; prolif-adjusted CDKN2A→KIF2C β=−0.098, p=1.1e-10) and not the
  CDK4/6–RB axis (co-dep with KIF18B r=0.55; CDK4/CDK6/RB1 r≈0). Obvious mechanism excluded; association real but unexplained.
- **KIF18A contrast.** Uniformly essential regardless of CDKN2A status (d=−0.09, n.s.) — no genotype selectivity; its named
  nomination (0.99→0.00 masked) was name recognition. The two CIN kinesins are not co-equal.
- **Novelty (5-angle literature prior-art sweep).** Four layers: (1) KIF2C as glioma target — NOT novel (prognostic since 2011;
  2023 inhibitor screen); (2) KIF2C in chromosomal instability — NOT novel (Compton/Bakhoum MCAK-CIN, 2007–2010); (3) 9p21-loss
  SL as a GBM strategy — NOT novel (active field; known nodes PRMT5 / CDK4-6 / TYMS-Chk1); (4) **KIF2C specifically as the
  CDKN2A/9p21 SL partner — NOVEL, no prior report found.** Effect-size context: KIF2C (d=−0.46 to −0.70) sits below pipeline
  anchors WDR77←9p21 (d=−0.87) and E2F3←RB1 (d=−1.53). Framed as a biomarker-stratified hypothesis for isogenic validation, not
  a demonstrated SL. Artifacts: `kif2c_cdkn2a_novelty_ledger.csv`, `kif2c_cdkn2a_effectsize_context.csv` (from the candidate-biology sub-analysis).

## 12. Key decisions & provenance
- **Novelty ruling (original, upheld):** an approved GBM drug does NOT down-weight a candidate; novelty is a separate
  reported axis. The user reaffirmed this at Step 10 by choosing P2-soft (no gating). *Provenance correction: an earlier
  agent paraphrase attributed a "down-weight rediscoveries later" quote to the user; a full-transcript recall check found no
  such statement — the recorded ruling is no-down-weight. The paraphrase is retracted and not used.*
- **Single model, no consensus** (user framing): all LLM steps are one model.
- **B's contribution:** novelty-label bug-fix (verified) + the variance-influence insight (nominal weight ≠ ranking influence:
  dependency 0.24 weight but ~6% influence because all candidates are pre-selected dependencies; selectivity 33%, tractability 21%).
- **PDCD5 directionality stress-test** carried from Step 4 through Step 10: ranked #127 on its own; Opus adjudication flagged it as a
  pro-apoptotic tumor suppressor needing reactivation, not inhibition. Documented as a design validation, not a nomination.
