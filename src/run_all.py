#!/usr/bin/env python
"""
run_all.py — run the full GBM target-discovery pipeline end to end.

Usage
-----
    python src/run_all.py                 # all steps, cached mode (offline)
    python src/run_all.py --from 08       # resume from a given step
    python src/run_all.py --only 09 10    # run just these steps
    python src/run_all.py --list          # list the steps and exit
    GBM_LIVE=1 python src/run_all.py      # recompute connector/LLM steps live

Modes (see src/common.py)
-------------------------
* cached (default): every step reproduces the committed results from data/cache/
  and results/ without network access or API keys. Steps 1/3/4/5/6 still need the
  large raw inputs (run step00 first) to recompute; the LLM steps (9/10/11a) load
  their committed run tables and need nothing external.
* live (GBM_LIVE=1): step00 re-downloads raw inputs and the GTEx connector; steps
  9/10/11a re-query the Anthropic API (needs the 'anthropic' package and
  ANTHROPIC_API_KEY). LLM steps are stochastic and will not reproduce the
  committed values bit-for-bit.

Iterative note: the final ranking depends on the gene-masked ablation (step11a)
and the expression-SL scan (step03). In cached mode any order works because each
step loads committed inputs. For a full LIVE recompute, run step11a before the
final ranking assembly in step10 (this script already orders 11a before a second
pass of 10 when --live-two-pass is given).
"""
import argparse
import os
import runpy
import sys
import time

SRC_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SRC_DIR)

# Ordered pipeline. Keys are the step tags used by --from / --only.
STEPS = [
    ("00", "step00_fetch_data",             "Data acquisition"),
    ("01", "step01_dependency",             "Glioma-selective dependency (DepMap)"),
    ("03", "step03_synthetic_lethality",    "Synthetic lethality"),
    ("04", "step04_tumor_vs_tme",           "Tumor-vs-TME deconvolution"),
    ("05", "step05_therapeutic_window",     "Therapeutic window / safety"),
    ("06", "step06_tumor_enrichment",       "Tumor expression enrichment (TCGA-GBM)"),
    ("07", "step07_tractability",           "Tractability + druggability + novelty"),
    ("08", "step08_feature_matrix",         "Feature matrix"),
    ("09", "step09_weight_selection",       "LLM weight selection"),
    ("11a", "step11a_gene_masked_ablation", "Gene-masked nomination ablation"),
    ("10", "step10_nomination",             "LLM nomination + final ranking"),
    ("11", "step11_ablations",              "Ablations + validation"),
    ("11b", "step11b_kif2c",                "KIF2C 9p21/CDKN2A SL deep-dive"),
]
# Note: 11a is ordered before 10 so a live recompute has the masked-ablation
# columns available when step10 assembles the final ranking.


def main():
    ap = argparse.ArgumentParser(description="Run the GBM target-discovery pipeline.")
    ap.add_argument("--from", dest="start", metavar="STEP", help="resume from this step tag (e.g. 08)")
    ap.add_argument("--only", nargs="+", metavar="STEP", help="run only these step tags")
    ap.add_argument("--list", action="store_true", help="list steps and exit")
    args = ap.parse_args()

    if args.list:
        live = os.environ.get("GBM_LIVE", "0") not in ("0", "", "false", "False")
        print(f"Mode: {'LIVE' if live else 'cached (offline)'}\n")
        for tag, mod, desc in STEPS:
            print(f"  [{tag:>3}] {mod:32s} {desc}")
        return

    tags = [s[0] for s in STEPS]
    if args.only:
        selected = [s for s in STEPS if s[0] in set(args.only)]
    elif args.start:
        if args.start not in tags:
            ap.error(f"unknown step {args.start!r}; choose from {tags}")
        i = tags.index(args.start)
        selected = STEPS[i:]
    else:
        selected = STEPS

    live = os.environ.get("GBM_LIVE", "0") not in ("0", "", "false", "False")
    print(f"Running {len(selected)} step(s) in {'LIVE' if live else 'cached (offline)'} mode.\n")

    failures = []
    for tag, mod, desc in selected:
        path = os.path.join(SRC_DIR, mod + ".py")
        print(f"\n{'=' * 72}\n[{tag}] {desc}\n{'=' * 72}")
        t0 = time.time()
        try:
            runpy.run_path(path, run_name="__main__")
            print(f"  [{tag}] done in {time.time() - t0:.1f}s")
        except Exception as exc:  # noqa: BLE001 - report and continue
            print(f"  [{tag}] FAILED: {type(exc).__name__}: {exc}")
            failures.append((tag, desc, f"{type(exc).__name__}: {exc}"))

    print(f"\n{'=' * 72}\nPipeline finished. {len(selected) - len(failures)}/{len(selected)} steps OK.")
    if failures:
        print("Failed steps:")
        for tag, desc, err in failures:
            print(f"  [{tag}] {desc}: {err}")
        sys.exit(1)


if __name__ == "__main__":
    main()
