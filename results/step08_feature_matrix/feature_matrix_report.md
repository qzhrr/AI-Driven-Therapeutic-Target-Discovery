# Step 8 — Feature matrix assembly + normalization

## Objective
Merge all six evidence layers into one row-per-candidate matrix (283 genes), every scored feature oriented
**higher = more favorable for a GBM small-molecule inhibition target**, normalized to [0,1], ready for LLM
category weighting (Step 9).

## Merge
Each layer's source file is authoritative for its native features (Step 5 re-pulled a few earlier columns;
deduped to original source). Joined on gene symbol → 283/283 rows, 0 NaN in scored features after imputation.

## 12 scored features → 7 categories
| Category | Features | Transform |
|---|---|---|
| Dependency | f_dep_selectivity (−Cohen's d, rank), f_dep_breadth (glioma_dep_fraction) | rank-pct / [0,1] |
| Selectivity | f_pan_ess_penalty (1 − pan_ess_fraction) | inverted |
| Synthetic-lethality (bonus) | f_sl_strength (sl effect), f_sl_breadth (n contexts) | min-max, 0-preserving |
| Tumor specificity | f_sc_specificity (rank), f_sc_tme_penalty (TME-dominated→0) | rank / binary |
| Tumor enrichment | f_tumor_enrichment (log2FC, rank), f_tumor_prevalence (pct over normal) | rank / [0,1] |
| Therapeutic window | f_safety (safety_composite) | [0,1] |
| Tractability | f_tractability | [0,1] |
| **Novelty (separate axis)** | f_novelty | [0,1] |

Normalization: rank-percentile for unbounded continuous features; already-[0,1] features kept as-is; SL
min-max (preserves the 0 for the 216 candidates with no SL). Penalties inverted so higher = better.

## Missing data
5 genes absent from GBmap (CSH2, KRTAP21-1, NOPCHAP1, POTEI, TUBB) → SC features imputed to neutral 0.5 +
`sc_measured=False`. All other features complete.

## Non-scored flags carried for the LLM (context, not weight)
has_inhibitor_moa, is_SL_partner_of_lost_TSG, lost_TSG_partner, sl_best_context, archetype, sc_compartment,
sc_specificity_tier, is_common_essential, uniprot_kinase/enzyme/membrane, uniprot_loc, chembl_action_types,
min_subgroup_n, ot_sm_approved_drug.

## Diagnostics (equal-weight baseline composite — NOT the final ranking)
- Range 0.207–0.746. Top-20: KIF2C, PDGFRA, C3orf38, CDK2, ITGB5, MET, JUN, KIF18B, WDR77, KIF18A...
- **Positive controls land well:** PDGFRA #2, WDR77 #9 (PRMT5 axis), MET #6, VRK1 #19, CDC25B #18.
- **E2F3 #14** — lifted by f_sl_strength=1.0 (strongest RB1-loss SL). SL bonus works.
- **PTK2 #100** — correctly demoted despite tractability=1.0, by poor safety (high normal-brain). Safety axis works.
- **PDCD5 #149 (middle):** the directionality stress-test. Pulled down by tractability=0 + poor safety, but NOT
  fully caught by numbers alone. `has_inhibitor_moa=False` is the flag the LLM must use in Step 10. Two-layer design confirmed.

## Feature independence (Spearman)
- sl_strength ~ sl_breadth ρ=0.98 (redundant → absorbed within SL category)
- dep_breadth ~ pan_ess_penalty ρ=−0.96 (intended tension, different categories)
- tumor_enrichment ~ tumor_prevalence ρ=0.88 (both TCGA → within enrichment category)
- tractability ~ novelty ρ=−0.53 (built-in tension, kept as separate axes)
Category-level weighting absorbs redundancy within categories and preserves tensions across them.

## Deliverables
- `feature_matrix.parquet` — 283 × 63 (raw + normalized + categories + baseline) — the master matrix
- `feature_matrix_scored.csv` — 12 features + 7 categories + baseline (readable)
- `feature_matrix.png` — feature correlation heatmap + baseline composite distribution
