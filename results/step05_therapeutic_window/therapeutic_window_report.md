# Step 5 — Therapeutic Window: normal-tissue + normal-brain safety filter

## Purpose
Filter/score the 283 glioma-selective dependencies by therapeutic window: a good small-molecule target is
expressed/required in tumor but LOW in normal tissues — especially normal brain, since GBM is intracranial and
on-target CNS toxicity is the dominant GBM-specific liability.

## Data sources (all 283 candidates covered)
1. **GTEx v8 body-wide** (54 tissues, `expression` connector, GENCODE v26 TPM). 274 candidates resolved directly;
   9 recovered via legacy symbols (POLR1H→ZNRD1, EPRS1→EPRS, NARS1→NARS, TARS1→TARS, H3C8→HIST1H3H, ODR4→C1orf27,
   POLR1G→CD3EAP, HJV→HFE2, NOPCHAP1→C12orf45). **283/283 covered.**
2. **GTEx brain** — 13 brain subregions (bulk depth).
3. **GSE97930 (Lake et al. 2018)** — dedicated healthy-brain snDrop-seq, 3 regions (frontal/visual cortex, cerebellum),
   ~35k cells, cell-type resolved. **Normal-neuron expression** (excitatory + inhibitory) is the key CNS-toxicity axis.
   265/283 measured. Cell-type assignment validated vs canonical markers (SNAP25/RBFOX3→neuron, AQP4/GFAP→astrocyte,
   MBP/PLP1→oligodendrocyte, CSF1R→microglia, CLDN5→endothelial — all correct).
4. **GBmap normal compartments** — oligodendrocyte/astrocyte/radial-glial/neuron/OPC from the Step-4 tumor atlas
   (tumor-context normal proxy).

## Scoring approach — GRADED, not hard-gated
A naive absolute normal-brain cutoff flags most candidates — a combined CNS flag (GTEx brain >10 TPM OR normal-neuron
in top quartile) hits 234/283, and the pure >10 TPM cutoff alone hits 229/283 (~81%) — GBM's brain location means almost
everything is "expressed in brain," so a hard gate nukes the list. Instead, each safety axis is **rank-normalized to
[0,1]** (1 = safest = lowest normal expression among candidates), and combined into `safety_composite`. Hard flags are
reserved for **top-decile risk only** (`cns_high_risk` n=47, `vital_high_risk` n=29) as LLM-visible annotations, not gates.

## Features added
- Normal-tissue: `gtex_max_tpm`, `gtex_max_tissue`, `gtex_vital_organ_max`
- Normal-brain: `gtex_brain_max` (bulk), `lake_neuron_max` / `lake_brain_max` (dedicated snRNA), `gbmap_normal_brain_max`
- Graded safety: `safety_vital_organ`, `safety_brain_bulk`, `safety_brain_neuron`, `safety_body_wide`, `safety_composite`
- Risk flags: `cns_high_risk`, `vital_high_risk`

## Results
- **Safest windows (tumor-expressed, low normal): UBE2C, ZNF593, LSM6, ACTL6A, PDCD2, RRP15** — cell-cycle/proliferation and
  small nuclear/processing genes with low normal-brain expression.
- **Riskiest (no window): ALDOA, RPS15, PFDN5, UBC, RPS11, RPS25** — ubiquitous ribosomal/housekeeping genes
  (RPS/RPL, ALDOA, UBC), extremely high in brain AND vital organs. Correctly down-weighted.
- Positive-control dependencies MET and VRK1 sit in the low-normal-brain region (good windows), consistent with
  their being clean glioma-selective dependencies.

## Reviewer carried-forward items (folded in here)
- **`min_subgroup_n`** (fragility): number of glioma lines supporting each dependency. Range 15–72; **zero candidates
  rest on <6 lines** (the ≥20% subset rule guarantees ≥~15 lines) — the reviewer's small-n concern does not materialize.
- **`sc_specificity_tier`**: confirmed-malignant-enriched (4) / not-TME-dominated-weak (231) / TME-dominated (43) /
  not-measured (5). Prevents treating weak positive specificity as strong evidence.
- **SC feature renamed** `sc_compartment` = tumor-cell-expressed (240) | TME-dominated (43) — an expression-compartment
  label, NOT a functional/oncogenic call (the PDCD5 lesson).
- **Directionality-unresolved list** seeded with PDCD5 (`directionality_unresolved.csv`), to be resolved at Steps 6+10.

## Deliverables
- `therapeutic_window_features.csv` — 283 candidates × 23 features (feeds the matrix)
- `gtex_safety_features.csv`, `lake2018_brain_features.csv` — per-source intermediates
- `directionality_unresolved.csv` — PDCD5 + future wrong-direction cases
- `therapeutic_window.png` — tumor-vs-normal-brain scatter + composite-safety distribution
