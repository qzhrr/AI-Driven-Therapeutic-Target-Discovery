# Gene-masked nomination ablation (Step 11 add-on)

**Question:** Does the LLM nomination layer nominate on the *evidence the pipeline provides*, or on
*gene-name recognition* (training-corpus familiarity)? This is the decisive test of the novelty thesis —
a "novel target discovery" pipeline is only credible if its nominations survive removal of the names.

**Design.** Two arms, identical harness (same model `claude-opus-4-8`, same system/tool/prompt, 100 runs
each), differing in ONE variable: in the **named** arm each candidate carries its real gene symbol; in the
**masked** arm the symbol is replaced by an opaque code (e.g. `CX657`) while every score, category value,
and pipeline-derived flag (SL partner, protein class, novelty tier, CNS/vital risk) is byte-identical.
A gene whose nomination frequency **holds** under masking is carried by the evidence; one that **drops** was
riding on its name; one that **rises** was being held back by its name.

Valid runs: named 100/100, masked 87/100 (13 masked runs errored; frequencies normalized per valid runs).
Zero hallucinated identifiers in either arm.

## Headline finding — the name matters more than the evidence

Across all nominated candidates, **masked frequency tracks the pipeline's own evidence far better than named
frequency does:**

| frequency correlates with… | NAMED ρ | MASKED ρ |
|---|---|---|
| selectivity (the top evidence category) | **0.12** | **0.57** |
| composite score (all evidence) | 0.47 | 0.72 |

named-vs-masked frequency correlate only ρ=0.32 — i.e. **the two arms disagree more than they agree.** When
you take the names away, the model's nominations snap back toward the quantitative evidence. That means in the
named arm, gene identity is *decoupling* nominations from the evidence — the definition of familiarity bias.

The fresh named arm reproduced the original Step-10 consensus core **exactly** (same 10 genes at ≥0.85),
confirming the ablation harness is faithful and the effect below is not harness drift.

## Three classes of behavior

**Evidence-carried (robust) — nominate on merit, name-independent (Δ≈0):**
VRK1 (1.00→0.99), MET (1.00→1.00), PDGFRA (1.00→1.00), CDK2 (1.00→0.91), ELAVL1 (0.98→0.92),
EGLN1 (0.95→0.99), plus E2F3, KIF18B, FGFR1. **These are the trustworthy nominations** — the model picks
them whether or not it knows what they are. Note this includes the approved-drug controls (MET, PDGFRA):
their nomination is evidence-driven, not just name-recognition.

**Name-INFLATED (famous → drops when masked) — the familiarity bias, made explicit:**
| gene | named | masked | Δ | why the drop is diagnostic |
|---|---|---|---|---|
| **KIF18A** | 1.00 | **0.00** | −1.00 | selectivity 0.11 — nominated on the *name* (clinical-stage inhibitors), collapses without it |
| **CDK7** | 1.00 | 0.03 | −0.97 | selectivity 0.00 — pure name recognition; no pipeline evidence supports it |
| PTPN11/SHP2 | 0.60 | 0.00 | −0.60 | famous allosteric target; CNS-risk + common-essential; gone when anonymized |
| KIF2C | 0.86 | 0.29 | −0.57 | partly name-driven, but retains some masked signal (selectivity 0.82) |
| TEAD1 | 1.00 | 0.47 | −0.53 | "genuinely druggable Hippo effector" — half the enthusiasm was the name |
| WDR77 | 0.62 | 0.22 | −0.40 | PRMT5-axis fame; selectivity 0.01 |

**Name-SUPPRESSED (rises when masked) — a prior working *against* the evidence:**
JUN (0.03→1.00), AMY2A (0.00→0.82), FERMT2 (0.00→0.53), C3orf38 (0.14→0.54). The model refused JUN when named
("undruggable AP-1 TF") but nominated it ~unanimously on its (strong: selec 0.83, dep 0.77, tract 1.00) numbers
alone. This is the mirror image of the bias — a *name-based veto* overriding favorable evidence. Whether that
veto is correct (AP-1 really is hard) or an over-correction is itself a judgment the numbers don't make.

## What this does to the four novel headline nominations

| gene | named | masked | Δ | verdict |
|---|---|---|---|---|
| **VRK1** | 1.00 | 0.99 | −0.01 | **robust — evidence-carried. The strongest novel pick survives cleanly.** |
| **ELAVL1** | 0.98 | 0.92 | −0.06 | **robust — evidence-carried.** |
| **KIF2C** | 0.86 | 0.29 | −0.57 | partly name-driven; retains selectivity signal but demote from headline. |
| **KIF18A** | 1.00 | 0.00 | −1.00 | **name-inflated — does NOT survive masking. Remove from the novel headline.** |

This confirms the concern raised in the results review from the evidence side: **VRK1 and ELAVL1 are the
defensible novel nominations; KIF2C is borderline; KIF18A was carried by name recognition, not by the
pipeline.** The two kinesins are not co-equal discoveries.

## Interpretation

1. **The LLM layer adds genuine value AND genuine bias — both are now measured.** The "promotions" the final
   report credited to expert reasoning (CDK7, TEAD1, PTPN11) are substantially name-recognition: they evaporate
   when anonymized. The report's charitable framing ("mechanistically strong targets the additive score buried")
   is only partly supported — for CDK7 (masked 0.03) it is not supported at all.
2. **The defensible claim is narrower and cleaner:** report VRK1, ELAVL1 (and the recovered controls
   PDGFRA/MET) as evidence-carried nominations; present KIF2C with the caveat; drop KIF18A from the headline
   or explicitly label it name-dependent.
3. **The masked ranking is arguably the more honest primary output** — it is the nomination the model makes
   from your evidence without importing outside priors. Consider reporting masked frequency as the headline
   and named as the "with-prior-knowledge" overlay, so the reader sees both.

## Deliverables
- `gene_masked_ablation.csv` — per-gene named vs masked frequency, Δ, effect class, evidence scores
- `gene_masked_ablation.png` — named-vs-masked scatter + Δ-frequency bars
- `ablation_named_raw.csv`, `ablation_masked_raw.csv` — all nominations, both arms

*Run 2026-07-13 15:10 UTC. Model claude-opus-4-8, 100 runs/arm, temperature default (deprecated on Opus).*
