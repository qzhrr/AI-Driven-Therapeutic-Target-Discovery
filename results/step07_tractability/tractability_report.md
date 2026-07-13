# Step 7 — Tractability + druggability annotation

## Purpose
Add the "is this actually druggable by a small molecule?" evidence layer, and — per the agreed design —
split it into **two independent features**: **tractability** (can a small molecule engage the target?) and
**novelty** (is it un-drugged?). Approved-drug status does NOT down-weight a candidate; it only lowers novelty.

## Four sources (all four in the therapeutic-window DB list)
1. **Open Targets** GraphQL — small-molecule (SM) tractability buckets (approved drug, advanced clinical,
   structure-with-ligand, high-quality pocket/ligand, druggable family) + drug/clinical-candidate count.
   All 283 candidates, 0 errors.
2. **ChEMBL** — target → mechanism (`action_type`: INHIBITOR/ANTAGONIST/AGONIST…) + genuine potent-compound
   count (IC50/Ki/Kd/EC50, pChEMBL≥6). Run on 61 drug/ligand-bearing genes + PDCD5.
3. **DGIdb** v5 GraphQL — drug-gene interactions (n drugs, n approved, interaction types). 164 genes returned,
   69 with ≥1 drug.
4. **UniProt** — protein class (kinase/enzyme), subcellular localization, membrane flag. 61 genes: 7 kinases,
   25 enzymes, 14 membrane.

## Feature definitions
- `tractability_score` (0–1, graded ladder): HQ pocket/ligand +0.35, structure-with-ligand +0.25,
  druggable family +0.20, ≥10 genuine potent compounds +0.20.
- `novelty_score`: 1.0 = undrugged, 0.5 = clinical/candidate, 0.0 = approved SM drug.
- `has_inhibitor_moa`: known INHIBITOR/ANTAGONIST/blocker chemistry (directionality signal).

## Data-quality fix (important)
ChEMBL's raw compound count inflated screening-only targets: **LMNA showed 6,192 "potent" compounds**, all
high-throughput screening `Potency` hits, not target-based IC50/Ki. Re-counting on genuine dose-response types
(IC50/Ki/Kd/EC50) collapsed LMNA 6,192→1 while real drug targets were unchanged (MET 6,729, FGFR1 5,414,
CDK2 2,409). The corrected count is what feeds `tractability_score`.

## Results
- **57 candidates have an approved SM drug** (novelty 0), 2 in clinical, 224 undrugged.
- **20 candidates are TRACTABLE + NOVEL** (tractability≥0.5, undrugged) — the sweet spot:
  CDC25B, VRK1, TARS1, FARSA, AK6, HSPA9, ALDOA, LONP1, OSBP, WDR77, ACTR3, GPX4, SNUPN, RABGGTB, CYFIP1, KIF18A, KIF2C, UBC, NCKAP1, ELAVL1.
- Standouts: **CDC25B** (phosphatase, 80 genuine compounds, druggable pocket), **VRK1** (kinase, structure+family,
  no compounds yet — genuinely novel), **WDR77** (PRMT5-axis SL gene), **GPX4** (ferroptosis), **KIF2C/KIF18A** (kinesins),
  **TARS1**.

## PDCD5 directionality test (first real check)
PDCD5 has a ChEMBL target entry (CHEMBL6066424) but **zero action types, zero mechanism, zero potent compounds** —
no small-molecule inhibitor has ever been developed. UniProt: nuclear/cytosolic, non-enzymatic. This is the empirical
fingerprint of a target you would want to *reactivate*, not inhibit: `has_inhibitor_moa=False`, `tractability_score`
low. The two-feature design + directionality flag catch what the pure expression/dependency numbers could not — exactly
the PDCD5 stress-test we set up in Step 4.

## Deliverables
- `tractability_features.csv` — 283 × 18 (feeds the matrix)
- `tractability_novelty.png` — (A) tractability×novelty axes, (B) 20 tractable+novel candidates, (C) chemical-matter landscape
