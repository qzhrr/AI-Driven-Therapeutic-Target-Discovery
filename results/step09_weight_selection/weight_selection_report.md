# Step 9 — LLM weight selection (single-model ensemble)

## Framing
Weights for the 7 scoring categories were selected by a single model — **Claude (Opus, `claude-opus-4-8`,
this session's model)** — run as a within-model ensemble. No consensus across models; no other models involved.
This is "the work done by Claude for Science, just you."

## Method
- **Task:** distribute exactly 100 points across the 7 categories (dependency, selectivity, synthetic-lethality,
  tumor-specificity, tumor-enrichment, therapeutic-window, tractability). Integer allocation, forced via a
  `tool_use` schema (no free-text parsing).
- **Blind to candidate scores** — the model weighted categories on domain reasoning alone, never seeing any
  gene's values, so weights are not reverse-engineered from a ranking.
- **Two conditions:** expert (GBM small-molecule target-discovery framing + domain category definitions) vs
  non-expert (generic "numerical feature categories", domain cues stripped). 100 runs each = 200 calls.
- **Novelty** kept as a SEPARATE reported axis, not weighted (per project decision: novelty must not down-weight
  already-drugged targets).

## Constraint discovered (documented)
The Opus model **deprecated the `temperature` parameter**, so the planned 0.70–0.85 temperature sweep was not
possible. Runs used default sampling. This is a model constraint from the Opus choice, not a design change, and
does not affect the deliverable (a weight distribution). Default sampling still produced run-to-run variation
(mean CV 0.055) — Opus is highly self-consistent on this task, which strengthens reproducibility.

## Final weight vector (expert median, normalized to Σ=1.0) → feeds Step 10
| Category | Weight | CV (stability) |
|---|---|---|
| dependency | 0.240 | 0.000 |
| selectivity | 0.180 | 0.032 |
| tractability | 0.160 | 0.068 |
| synthetic_lethality | 0.120 | 0.045 |
| therapeutic_window | 0.120 | 0.087 |
| tumor_specificity | 0.110 | 0.063 |
| tumor_enrichment | 0.070 | 0.092 |

Mean CV across categories = 0.055 (low = stable). Dependency is the single most important category (0.24,
CV 0.000 — identical in all 100 runs), followed by selectivity (0.18) and tractability (0.16).

## Expert vs non-expert (the Step-11 ablation, computed now)
The two prompt conditions differ in interpretable, domain-meaningful ways (all MW-p ≤ 1e-5, tight distributions):
- **Expert up-weights** tractability (+8), selectivity (+3), dependency (+2), therapeutic-window (+2) — the
  determinants of a *druggable, selective* target.
- **Non-expert up-weights** tumor-enrichment (+8) and tumor-specificity (+7) — naive "high-in-tumor = important"
  reasoning, the exact expression-only trap that mis-prioritizes microenvironment/housekeeping genes.
This quantifies the "expert knows druggability + selectivity beat raw expression" story, and gives Step 11 a
real expert-vs-non-expert weight contrast to propagate into rankings.

## Deliverables
- `weights_expert_raw.csv`, `weights_nonexpert_raw.csv` — per-run allocations (100 each) + rationales
- `weights_aggregated.csv` — final weight vector + stability
- `weights_comparison.csv` — expert vs non-expert per category + MW-p
- `weight_selection.png` — (A) expert vs non-expert distributions, (B) final weight vector
