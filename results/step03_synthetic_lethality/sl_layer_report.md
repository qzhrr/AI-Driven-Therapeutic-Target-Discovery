# Step 3 — Synthetic-Lethality / Genotype-Stratified Dependency Layer

**Data:** DepMap 24Q4 — CRISPR (Chronos) gene effect + OmicsAbsoluteCNGene + OmicsSomaticMutationsMatrixDamaging + OmicsSomaticMutationsMatrixHotspot.
**Scope:** 283 glioma-selective candidates × 9 GBM-relevant lesion contexts, 1,178 cell lines (72 glioma).

## Method (pan-cancer detection, glioma-anchored interpretation)

SL biology is largely lineage-agnostic, so associations are **detected pan-cancer** (full power, ~1,100 lines),
but only **count** when they also (a) hold in the same direction within glioma and (b) target a lesion prevalent
in glioma. Pan-cancer is the detection lens; glioma is the frame of meaning.

Per candidate × lesion: gene-effect altered vs wild-type, Cohen's d + Mann-Whitney U, BH-FDR over 2,547 tests.
SL-like = candidate **more** dependent in altered lines (d < −0.3, FDR < 0.05).

**Glioma-prevalence floor = 8%** (lowered from an initial 15% after review, to admit PDGFRA amplification — a
bona fide GBM oncogene context present in ~8.5% of glioma lines. Trade-off: the glioma-consistency check on an
8.5%-prevalence lesion rests on ~6 glioma lines, so those hits are direction-consistent but individually
lower-powered within glioma; they were detected pan-cancer where power is ample.)

**Lesion prevalence in glioma:** CDKN2A del 59.6%, CDKN2B del 58.5%, TP53 mut 67.0%, PTEN loss 48.9%,
MTAP del 38.3%, EGFR amp 23.4%, NF1 loss 17.0%, RB1 loss 16.0%, PDGFRA amp 8.5%. CN thresholds: homozygous
deletion CN < 0.5, amplification CN ≥ 6 (absolute copy number).

## Results

| Filter stage | Associations |
|---|---|
| Pan-cancer SL-like (d<−0.3, FDR<0.05) | 109 |
| + glioma-direction consistent | 91 |
| **+ glioma-relevant lesion, ≥8% (anchored)** | **91** |

**67 of 283 candidates** carry ≥1 anchored SL association; **48** are flagged as the synthetic-lethal
partner of a recurrently-lost GBM tumor suppressor.

**Anchored hits by lesion:** PTEN_loss 18, CDKN2A_del 14, CDKN2B_del 13, MTAP_del 12, PDGFRA_amp 12,
RB1_loss 8, EGFR_amp 6, NF1_loss 4, TP53_mut 4.

### Canonical SL relationships recovered (method validation)

- **E2F3 ← RB1 loss** (d = −1.53, strongest). RB1 loss de-represses E2F → E2F3 dependency. Textbook RB-pathway SL.
- **WDR77 ← 9p21 co-deletion (MTAP / CDKN2A / CDKN2B)** (d = −0.87). WDR77/MEP50 is the obligate cofactor of
  **PRMT5** — the MTAP-deletion→PRMT5-axis SL, the most clinically-advanced SL in oncology. Recovery validates the approach.
- **CDK2 / INTS8 ← RB1 loss** — RB/cell-cycle logic.
- **MLST8 / RPTOR / PTPMT1 / FOXK1 ← PTEN loss** — PTEN→PI3K/mTOR axis; MLST8 and RPTOR are mTOR-complex components.
- **PDGFRA-amp amplicon block** (RPTOR, RNF4, TUBB, MSTO1, PFDN4, VBP1, LSM6, PPP1CB, SLC35B2, EXT1, CHMP7,
  C3orf38) — dependencies enriched in PDGFRA-amplified lines; admitted at the 8% floor.
- **GRB2 ← EGFR amplification** — EGFR-signalling adaptor dependency.

## Interpretation & honest caveats

- **Not every candidate has an SL hit, and that is correct.** Lineage/adhesion-driven candidates (JUN, integrin
  cluster) have no genotype association; the SL layer is an *added bonus annotation* for the genotype-driven
  subset, not a filter all candidates pass.
- **Genotype-association SL ≠ true digenic SL.** DepMap is single-gene KO; this measures dependency enriched in
  altered lines — a strong proxy, not a double-knockout measurement.
- **9p21 co-deletion.** MTAP, CDKN2A, CDKN2B are adjacent on 9p21 and co-deleted as one event, so WDR77/PRMT5-axis
  hits appear against all three; the mechanistic driver is MTAP metabolic loss (PRMT5 axis).
- **PDGFRA-amp power.** Admitted deliberately at the 8% floor. The 12 genes are glioma-direction-consistent but,
  given ~6 amplified glioma lines, individually underpowered *within* glioma — they earned entry on pan-cancer
  detection + consistency, not on strong within-glioma statistics. Flagged for the LLM to weigh accordingly.
- **TSG classification.** Cross-checked vs Open Targets hallmarks: CDKN2A, NF1, RB1, TP53 = TSG; EGFR, PDGFRA =
  oncogene. Open Targets mislabels PTEN "oncogene" (known single-token artifact); PTEN handled by loss-of-function
  mechanism (deletion / damaging mutation), which is correct.

## Features added to the matrix

`sl_genotype_max_effect_size` (−Cohen's d of best anchored hit), `sl_best_context`, `sl_fdr`, `sl_n_contexts`,
`sl_all_contexts`, `is_SL_partner_of_lost_TSG` (flag), `lost_TSG_partner`, `sl_context_type`. Per the user's ruling,
SL-linkage enters scoring as a **bonus**, not a strong driver.

## Deliverables

- `sl_genotype_associations.csv` — all 2,547 candidate × lesion tests (anchoring at 8% floor)
- `sl_candidate_features.csv` — per-candidate SL features (feeds the feature matrix)
- `sl_anchored_associations.parquet` — 91 anchored hits (checkpoint)
- `sl_heatmap.png` — candidate × lesion SL-strength heatmap (top 30)

---

# SL layer v2 — expression-context detector (paralog buffering)

**Motivation (diagnosed in the B side-analysis).** The v1 module tests only 9 copy-number/mutation
lesions, so it is structurally blind to synthetic-lethal partners silenced *epigenetically or
transcriptionally* rather than deleted. The canonical example is **VRK1←VRK2**: VRK2 is rarely
deleted (v1 scored VRK1 `sl_effect=0.0`) but is promoter-methylated / expression-low across a large
fraction of gliomas — a validated paralog-buffering SL the genotype module could not see.

**Method.** Correlate each candidate's CRISPR dependency profile against genome-wide expression;
a paralog whose *low expression* predicts *stronger dependency* is flagged as a candidate SL partner
(`expression_SL_scan_candidates.csv`, 283×15; 44 significant hits, of which the mechanistically
interpretable class is paralog pairs — mostly buffering, plus one obligate co-complex co-dependency).

**Paralog verification (added here).** All 10 symbol-root paralog hits were checked against **Ensembl
Compara** (`ensembl_homology`, `homology_type=paralogues`, homo_sapiens). All 10 are genuine
within-species paralogs. Nine of the ten reach the significance threshold (ACTL6A←ACTL6B drops,
d=0.19, p=0.01), leaving **9 credible paralog SLs**, 6 of which v1 scored as exactly 0.0. Eight are
buffering pairs; **CHMP3←CHMP2A is included as a co-dependency SL** — CHMP3/CHMP2A are obligate
ESCRT-III co-complex members, so the relationship is complex-membership co-dependency rather than
paralog buffering, but it is a solid, significant dependency association and is retained on that
basis (per user decision):

| gene | partner | expr-SL \|Cohen's d\| | prior report |
|---|---|---|---|
| FAM50A | FAM50B | 3.14 | validated |
| CDS2 | CDS1 | 2.31 | validated |
| RPP25L | RPP25 | 2.16 | validated |
| EAF1 | EAF2 | 1.49 | validated |
| PTK2 | PTK2B | 1.16 | partial |
| ELMO2 | ELMO3 | 1.16 | validated |
| CHMP3 | CHMP2A | 0.83 | co-complex (ESCRT-III co-dependency, not buffering) |
| GOLT1B | GOLT1A | 0.70 | none (novel pair) |
| VRK1 | VRK2 | 0.70 | validated |

**Integration.** `sl_effect = max(genotype SL, expression-paralog SL)`; `sl_n_contexts += 1` when the
expression-SL is a new/stronger signal. Normalization: **winsorized at the strongest genotype SL
(E2F3, d=1.53)** before min-max, so `f_sl_strength` stays anchored to a biologically-benchmarked scale
and the FAM50A d=3.14 outlier does not compress genotype-SL genes (user's feature-score decision).

**Downstream impact.** Whole-ranking Spearman v1↔v2 = **0.993**. Top-100 membership changes by exactly
one gene (FAM50A in at #45, HAUS2 out). CHMP3 rises #255→#213 — deep in the list, so it does not touch
the reviewed pool. The headline is unaffected and slightly strengthened: **KIF2C stays #1**, and
**VRK1 rises #9→#7** now carrying independent expression-SL corroboration (VRK1←VRK2, validated) on top
of its dependency evidence. Steps 4–7 are computed independently of SL and are unchanged; Step-9 weights
were assigned blind and are unchanged; Step-11 ablation already ranked SL the least load-bearing
category (LOO ρ=0.981).

## v2 deliverables
- `expression_SL_scan_candidates.csv` — genome-wide expression-SL scan (from the B side-analysis)
- `sl_augmentation.png` — 3-panel: recovered paralog SLs, rank movement, overall stability
- `sl_v2_rank_delta.csv` — per-gene v1→canonical rank delta with SL provenance
- `feature_matrix.parquet` (v4), `final_candidate_ranking.csv` (v4) — canonical, winsorized SL-v2 (CHMP3 included)
