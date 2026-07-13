# Step 10 — LLM target nomination (single-model ensemble)

## Framing & posture
Nominations by a single model — **Claude (Opus, `claude-opus-4-8`, this session's model)** — no consensus across
models. **Posture P2-soft** (user decision): all 283 candidates in ONE ranked list by expert-weighted composite,
corrected 3-tier novelty as a visible tag; NO gating. Faithful to the original ruling that novelty is reported,
not down-weighted.

## Upstream correction adopted (from B's review, independently verified)
The Step 7/8 `novelty_score` mislabeled 57 genes as "drugged" via `ot_sm_approved_drug OR dgidb_n_approved_drugs>0`.
Only **21** carry the authoritative Open Targets approved-SM-drug flag; the other 36 were flagged by DGIdb alone —
all 36 simply had *any* DGIdb interaction record (promiscuous approved drugs create spurious "approved" links to
undrugged genes: JUN, PTPN11, TEAD1, VCP, RAC1, GRB2, PTK2). Replaced with a 3-tier axis on the authoritative
signal: **21 approved_sm_drug / 15 has_sm_chemistry / 247 novel**. This corrected the apparent rediscovery rate in
the top-20 from 9/20 to **4/20**. The correction does NOT change any composite score (novelty was always a separate
axis) — only tier labels and the presentation.

## Method
- Base: all 283 ranked by expert-weighted composite (Step-9 weights). LLM reviewed the **top 100** each run.
- Per candidate the model saw: 7 category sub-scores + composite + context flags (novelty_tier, has_inhibitor_moa,
  SL partner, sc_compartment, CNS-safety, protein class). Structured `tool_use` output.
- **100 runs**, default sampling (temperature deprecated on Opus). Each nomination carries an intended **modality**
  (inhibit/degrade/reactivate/unknown) + confidence → directionality tested on everything nominated.
- 1,497 total nominations, 0 hallucinated genes, 0 errors.

## Result — consensus core (nominated in ≥85 of 100 runs)
CDK2, CDK7, EGLN1, ELAVL1, KIF18A, KIF2C, MET, PDGFRA, TEAD1, VRK1

- **Validated GBM oncogene controls recovered:** PDGFRA (#2, freq 1.00), MET (#4, 1.00), FGFR1 (#11, 0.50) — the
  pipeline recovers known GBM biology.
- **Genuinely NOVEL consensus nominations:** **ELAVL1, KIF18A, KIF2C, VRK1** (all freq ≥0.87) — the headline output.
  KIF2C, KIF18A (mitotic kinesins, CIN-selective), VRK1 (paralog-lethal kinase), ELAVL1/HuR (RNA-binding oncoprotein).
- has-chemistry consensus: CDK2, CDK7, TEAD1 (druggable, not yet approved in GBM).
- **Modality reasoning works:** E2F3 nominated "degrade" 33× (undruggable TF → degrader), WDR77 "degrade" 7×.
  1,435/1,497 inhibit, 51 degrade, 11 unknown.

## PDCD5 directionality stress-test — PASSED (the decisive validation)
PDCD5 ranked #127 (outside top-100) on its own — its near-zero tractability + poor safety correctly sink it. In a
dedicated adjudication, Opus independently identified it as a **pro-apoptotic tumor suppressor LOST in GBM**, flagged
the **directionality error** (a small-molecule inhibitor would suppress a tumor suppressor — opposite of therapeutic
benefit; it needs REACTIVATION, not inhibition), and diagnosed that the tumorSpec=1.00/dependency scores could
mislead. The two-layer design works: the composite alone could not catch this; the LLM layer did.
Adjudication excerpt:
> # PDCD5 Assessment for GBM Small-Molecule Inhibition
> ...small-molecule *inhibitor* would suppress a pro-apoptotic tumor suppressor — the opposite of a therapeutic effect.

## Deliverables
- `final_candidate_ranking.csv` — 283 genes: composite + all category scores + LLM nomination frequency/rank + modality + novelty tier + control class
- `nomination_frequency.csv` — per-gene frequency, mean rank, dominant modality
- `nomination_runs_raw.csv` — all 1,497 nominations (run, rank, modality, reason)
- `nomination_frequency.png` — frequency top-20 (tier-colored) + consensus-vs-composite
- `nomination_controls.png` — GBM oncogene controls recovered + novel nominations
