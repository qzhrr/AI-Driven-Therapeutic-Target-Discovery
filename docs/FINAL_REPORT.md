# GBM small-molecule target discovery — Final Report
*Cross-dataset integration + single-model LLM weighting & nomination. All work by Claude for Science.*
*v4 — incorporates the gene-masked nomination ablation (blinded test), the KIF2C mechanism deep-dive, the
KIF2C–CDKN2A/9p21 synthetic-lethality prior-art assessment, and an upgraded SL layer (expression-context
detector for paralog/co-complex synthetic lethalities; ranks reflect this SL-augmented composite).*

## Executive summary
Starting from **17,916 genes**, a genetic-dependency core followed by a therapeutic-window filter and single-model
LLM weighting/nomination produced **283 glioma-selective dependency candidates**, ranked by an expert-weighted
7-category composite and 100 runs of LLM nomination. The pipeline **recovers known GBM biology** (positive controls
PDGFRA/MET/VRK1/FGFR1 in the top 12 and PTPN11 correctly down-weighted as common-essential; canonical false-positive EGFR correctly excluded). A **gene-masked
ablation** (nominating on anonymized candidates) separated evidence-carried nominations from name-recognition bias,
narrowing the novel headline to **VRK1 and ELAVL1**.

**The single most distinctive lead is a genotype-stratified one: KIF2C is a synthetic-lethal dependency of 9p21.3
(CDKN2A/CDKN2B/MTAP) deletion** — the most common focal deletion in IDH-wt GBM (glioma d=−0.70, Mann-Whitney p=8.3e-3;
pan-cancer d=−0.46, raw p=2.8e-10, BH-FDR=7.9e-8). A five-angle literature sweep found **no prior report of the
KIF2C–CDKN2A/9p21 pairing**, even
though the 9p21-SL field is active (its known nodes are PRMT5, CDK4/6, TYMS — not KIF2C). This is a biomarker-defined
hypothesis requiring isogenic validation, not a demonstrated SL — but it is the pipeline's clearest genuinely-new idea.

## The funnel
17,916 genes → (glioma-selective dependency, DepMap) → **283 candidates** → (6 evidence layers, 12 features, 7 categories)
→ expert-weighted composite → (100 LLM nomination runs, named + masked) → **evidence-carried novel core**.

## Headline result — the defensible novel targets (survive gene masking)
A nomination is only credible for a *novel-target* pipeline if it survives removal of the gene name. Re-running the
100-run nomination on anonymized candidates (all scores byte-identical, symbols replaced by opaque codes) gives:

All frequencies below are from the **ablation's own two arms** (named vs masked, same 100/87-run design), so named + Δ = masked
is internally consistent. (The Step-10 nomination frequencies used elsewhere in this report — VRK1 0.99, ELAVL1 1.00,
KIF18A 0.99 — are a separate 100-run batch and differ by ≤0.02.)

| Gene | ablation named freq | **masked freq** | Δ | verdict |
|---|---|---|---|---|
| **VRK1** | 1.00 | **0.99** | −0.01 | **robust — evidence-carried; strongest novel pick** |
| **ELAVL1** (HuR) | 0.98 | **0.92** | −0.06 | **robust — evidence-carried** |
| KIF2C (MCAK) | 0.86 | 0.29 | −0.57 | borderline — partly name-driven, but retains a real 9p21-del SL signal (see below) |
| KIF18A | 1.00 | **0.00** | −1.00 | **name-inflated — removed from headline; nominated on name, not evidence (selectivity 0.11)** |

**The two novel targets that stand on the evidence alone are VRK1 and ELAVL1.** Two more novel genes are also
evidence-carried at lower frequency: **E2F3** (0.63→0.79, *rises* when masked) and **KIF18B** (0.57→0.63). Masked-arm
frequencies are normalized over 87 valid masked runs (the batch completed runs 0–86; runs 87–99 are absent as a contiguous
block, i.e. batch truncation rather than scattered per-run failures — so the surviving runs are the first 87, not a biased
subset). Verdicts are stable across the 87 runs.

## KIF2C — a 9p21/CDKN2A-deletion-selective dependency (the pipeline's most distinctive lead)
KIF2C is borderline on name-masking (0.87→0.29), so it does **not** earn a headline slot on nomination frequency. What
makes it worth carrying is orthogonal to the LLM: **KIF2C dependency is selective for loss of the 9p21.3 locus** —
the single most common focal deletion in IDH-wt GBM.

- **Genotype-selective dependency.** In the 72-line glioma cohort, CDKN2A-deleted lines are markedly more KIF2C-dependent
  than intact lines (mean gene-effect −0.503 vs −0.353; **Cohen's d=−0.70, Mann-Whitney p=8.3e-3**, n=44 del / 28 intact;
  45% of deleted lines strongly dependent vs 14% intact). The signal replicates pan-cancer as an anchored synthetic-lethal
  association (**d=−0.46, raw p=2.8e-10, BH-FDR=7.9e-8, n=298** CDKN2A-deleted lines).
  **Attribution caveat:** in this cohort CDKN2A, MTAP and CDKN2B deletions are effectively inseparable (CDKN2A_del vs MTAP_del
  r=0.97, one discordant line of 72), and MTAP-deletion actually tracks KIF2C dependency marginally *more* strongly (d=−0.74)
  than CDKN2A (d=−0.70). No single-gene attribution within the 9p21 block is statistically supported — the biomarker is
  **9p21-block deletion**, not CDKN2A specifically. Within the block the per-lesion effect sizes are **MTAP d=−0.74**,
  **CDKN2A d=−0.70**, **CDKN2B d=−0.37** — i.e. it tracks the whole **9p21 co-deletion block**, not one gene. In this
  cohort 61% of lines are CDKN2A-deleted and all 44 are IDH-wild-type, so the biomarker maps onto the dominant GBM subtype.
- **Not proliferation, not the CDK4/6–RB axis.** The dependency does not track mitotic index (glioma r=+0.21 n.s.) and the
  proliferation-adjusted CDKN2A→KIF2C effect is unchanged (β=−0.098, p=1.1e-10); KIF2C co-depends with KIF18B (r=0.55), not
  CDK4/CDK6/RB1 (r≈0). The obvious mechanism is ruled out — the association is real but mechanistically unexplained.
- **KIF18A is the opposite** — uniformly essential regardless of CDKN2A status (d=−0.09, n.s.); no genotype selectivity, and
  its named nomination (0.99→0.00 masked) was clinical-stage-inhibitor name recognition. The two CIN kinesins are not co-equal.

### How novel is the KIF2C–CDKN2A/9p21 synthetic lethality? (prior-art assessment)
The claim is **novel at the level that matters — the specific pairing — but must be stated narrowly.** A five-angle
literature sweep separates four layers:
1. *KIF2C as a glioma target* — **not novel** (prognostic marker since 2011; a 2023 structure-based screen already pursued KIF2C inhibitors for glioma).
2. *KIF2C in chromosomal instability* — **not novel** (MCAK suppresses CIN when present / promotes it when lost; Compton/Bakhoum, 2007–2010).
3. *9p21-loss synthetic lethality as a GBM strategy* — **not novel** (active field; but the known nodes are **PRMT5** (MTAP arm), **CDK4/6** (CDKN2A arm), **TYMS/Chk1** — KIF2C is not among them).
4. *KIF2C specifically as a 9p21-block-deletion SL partner in GBM* — **novel; no prior report found.** (Because CDKN2A/MTAP/CDKN2B are collinear here, the claim is made about the 9p21 block as a whole, not CDKN2A alone.) This conjunction is the pipeline's actual contribution: it surfaced a 9p21-loss vulnerability node the existing SL literature has not considered.

**Caveats that cap the claim (do not overstate):** this is a **computationally-derived, unvalidated genotype association**,
not a demonstrated synthetic lethality — single-gene-KO dependency enriched in deleted lines (44 vs 28; underpowered), with
no isogenic CDKN2A knockout/rescue done. The **effect is modest** next to the pipeline's own anchors (WDR77←9p21 d=−0.87;
benchmark E2F3←RB1 d=−1.53). And the mechanism is unexplained (the obvious CDK4/6–RB route was excluded). It should be
presented as **a biomarker-stratified hypothesis requiring isogenic validation**, not an established SL. (Contrast VRK1:
*de-risked-but-not-novel* — its VRK2-paralog SL is already published; KIF2C–CDKN2A is *novel-but-unvalidated*.)

## ELAVL1 novelty caveat (real-world chemistry)
ELAVL1/HuR is evidence-robust but **not a chemically untouched target**. HuR has an active molecular-glue/PROTAC degrader
literature: "Druglike Molecular Degraders of the Oncogenic RNA-Binding Protein HuR" (*JACS Au*, doi:10.1021/jacsau.5c00551)
and "Molecular glue degraders of HuR suppress BRAF-mutant colorectal cancer" (*Nature* 2026, doi:10.1038/s41586-026-10613-5).
The dominant tractable modality is therefore **targeted degradation**, not classical small-molecule inhibition, and the
disclosed disease context (BRAF-mutant CRC) is not GBM. So ELAVL1 is best framed as a *degrader* opportunity that still
needs a CNS/GBM-specific path, whereas **VRK1 is the more genuinely novel small-molecule (inhibitor-lane) pick.**

## Positive-control recovery (validation)
- **4/5 pre-declared controls in the top 12:** PDGFRA (#2), MET (#4), VRK1 (#7), FGFR1 (#12); the fifth, PTPN11 (#74), is
  flagged common-essential and correctly down-weighted by the pan-essentiality penalty — expected behavior for a
  selectivity-weighted composite, not a miss. Mann-Whitney p = 3.1e-5 (declared controls ranked better than the rest).
  Controls MET/PDGFRA are **evidence-carried under masking** (Δ≈0) — recovered on merit,
  not just name.
- **Negative-direction validation:** EGFR, CDK4, CDK6 correctly did NOT survive Step 1 (EGFR = CRISPR de-amplification
  artifact, rank 17,834; CDK4/6 common-essential, not glioma-selective).

## What the LLM layer adds — value AND measured bias
The masked ablation quantifies both sides:
- **Genuine value:** masked frequency tracks the pipeline's own evidence *better* than named frequency (selectivity
  correlation ρ 0.12 named → **0.57 masked**; composite ρ 0.47 → 0.72). Stripped of names, the model snaps back to the evidence.
- **Measured bias:** some "expert promotions" in the v1 report were substantially name recognition. **CDK7** (named 1.00
  → masked 0.03, selectivity 0.00) was *not* supported by pipeline evidence — the v1 "mechanistically strong target the
  composite buried" framing was wrong for CDK7. **TEAD1** (1.00→0.47) and **PTPN11/SHP2** (0.59→0.00) were half-to-mostly name.
- **Name-based veto (mirror image):** JUN was refused when named ("undruggable AP-1 TF", 0.03) but nominated ~unanimously
  when masked (1.00) on strong numbers (selec 0.83, dep 0.77, tract 1.00). Whether the veto is correct or an over-correction
  is a judgment the numbers do not make.

## Ablations & robustness
composite-only vs +LLM (ρ=0.49/0.55); expert vs non-expert weights (ρ=0.917 — expert favors druggable, non-expert favors
proliferation); weight jitter ±0.03 (ρ=0.975, 89% top-20 retention); leave-one-out (selectivity + tractability load-bearing);
gene-masked ablation (named vs masked; see above).

## Caveats
- **2D-culture dependency** (DepMap CRISPR); **cohort skew** (72 lines, mostly IDH-wt GBM); **snDrop-seq depth** (rank signal
  only); **BBB/CNS penetrance deferred** to chemistry; **single model** (no consensus).
- **Name-recognition bias is real and now measured** — treat named nomination frequency as "with-prior-knowledge" and masked
  frequency as the evidence-only readout. Report both.
- **Masking controls for name, not for score-profile identity.** The tractability and novelty features are built from
  real-world drug data (ChEMBL potent-compound counts, Open Targets approved-SM flags), so a maxed score vector still leaks
  "known druggable target" even with the symbol hidden (e.g. MET masked→1.00 has tractability=1.0). Masked frequency is
  therefore a **lower bound** on identity influence, not a pure evidence readout — though for the novel picks the leak is
  small (VRK1 tractability 0.80, ELAVL1 0.55, not maxed), so their robustness survives.
- **Nomination is scoped to the composite's own top-100.** Each nomination run reviews only the top 100 by expert composite,
  so 183/283 candidates can never be nominated and nomination frequency is a re-ranking *within* the composite's selection,
  not an independent signal over the full list. Weight selection (Step 9) was firewalled (weights assigned blind to scores),
  so this is not circular, but the "semi-independent layer" adds signal only within that gated set.

## Bottom line — the portfolio, by what kind of lead each is
Three leads survive scrutiny, and they are novel/de-risked in *different* ways — the distinction is the actual result:
- **KIF2C — the distinctive novel lead: a 9p21/CDKN2A-deletion synthetic lethality.** Genotype-stratified (biomarker =
  the commonest IDH-wt GBM deletion), and the specific KIF2C–9p21 pairing has no prior art. Novel-but-unvalidated: a
  computationally-derived, modest-effect, mechanistically-unexplained correlation that needs isogenic CDKN2A knockout/rescue
  before it is a demonstrated SL. This is the pipeline's clearest new hypothesis.
- **VRK1 — the de-risked lead (#7).** Evidence-carried under masking (Δ≈−0.01) and a positive control; its VRK2-paralog SL is
  already published — and the upgraded SL layer now recovers it *within the pipeline* (VRK1←VRK2, d=−0.70), lifting VRK1 #9→#7.
  Validated-but-not-novel, and the most genuinely novel *small-molecule inhibitor* pick of the set.
- **ELAVL1/HuR — evidence-carried but chemically occupied.** Robust under masking, but HuR has published molecular-glue/PROTAC
  degrader chemistry (JACS Au; Nature 2026) — a degrader opportunity, not an untouched target.

The framework's value is precisely this triage: dependency-first evidence + a blinded name-masking test + a prior-art
check together separate *novel-and-unvalidated* (KIF2C) from *validated-and-not-novel* (VRK1) from *chemically-occupied*
(ELAVL1) — a distinction raw nomination frequency alone cannot make.

## Deliverables
`final_candidate_ranking.csv` (283 genes: all scores + composite + named & masked LLM frequency + Δ + effect class +
modality + novelty tier). See `README_deliverables.md` for the full index, including the masked-ablation, KIF2C-mechanism,
and KIF2C–CDKN2A novelty-ledger artifacts.
