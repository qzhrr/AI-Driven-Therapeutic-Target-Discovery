"""
common.py — shared paths and I/O helpers for the GBM target-discovery pipeline.

All pipeline scripts import from here so that data locations are defined in exactly
one place and every path is repository-relative (no absolute paths anywhere).

Layout (relative to the repository root):
    data/raw/     large external inputs — git-ignored, fetched by step00_fetch_data.py
    data/cache/   small connector/LLM responses — committed, enable offline reproduction
    results/stepXX_name/   per-step outputs (tables, figures, reports)

Run modes
---------
The two pipeline stages that are not pure computation — connector queries
(``host.mcp``) and LLM calls (``host.llm``) — read from data/cache/ by default so
the whole pipeline is reproducible offline. Set the environment variable
``GBM_LIVE=1`` to force those stages to recompute against the live APIs instead
(requires network access, connector availability and an Anthropic API key; see
adapters.py). Stochastic LLM steps (9, 10, 11a) can only be reproduced *exactly*
from the committed cache — a live re-run draws new samples.
"""
from __future__ import annotations

import os
from pathlib import Path

# ---------------------------------------------------------------------------
# Repository-relative paths
# ---------------------------------------------------------------------------
# src/ -> repository root is one level up.
SRC_DIR = Path(__file__).resolve().parent
ROOT = SRC_DIR.parent

DATA = ROOT / "data"
DATA_RAW = DATA / "raw"
DATA_CACHE = DATA / "cache"
RESULTS = ROOT / "results"

# Per-step result directories (created on import so scripts can write freely).
STEP_DIRS = {
    "step01_dependency": RESULTS / "step01_dependency",
    "step03_synthetic_lethality": RESULTS / "step03_synthetic_lethality",
    "step04_tumor_vs_tme": RESULTS / "step04_tumor_vs_tme",
    "step05_therapeutic_window": RESULTS / "step05_therapeutic_window",
    "step06_tumor_enrichment": RESULTS / "step06_tumor_enrichment",
    "step07_tractability": RESULTS / "step07_tractability",
    "step08_feature_matrix": RESULTS / "step08_feature_matrix",
    "step09_weight_selection": RESULTS / "step09_weight_selection",
    "step10_nomination": RESULTS / "step10_nomination",
    "step11_ablations": RESULTS / "step11_ablations",
    "step11a_gene_masked_ablation": RESULTS / "step11a_gene_masked_ablation",
    "step11b_kif2c": RESULTS / "step11b_kif2c",
}


def ensure_dirs() -> None:
    """Create the data and results directories if they do not yet exist."""
    for d in (DATA_RAW, DATA_CACHE, RESULTS, *STEP_DIRS.values()):
        d.mkdir(parents=True, exist_ok=True)


ensure_dirs()

# ---------------------------------------------------------------------------
# Run-mode flag
# ---------------------------------------------------------------------------
LIVE = os.environ.get("GBM_LIVE", "0") not in ("0", "", "false", "False")


def use_cache(path: os.PathLike | str) -> bool:
    """
    True when a cached artifact should be used instead of recomputing.

    A cached result is used when it exists on disk and live mode is off. This is
    the guard placed in front of every ``host.mcp`` / ``host.llm`` block so the
    pipeline reproduces offline yet can be forced live with GBM_LIVE=1.
    """
    return (not LIVE) and Path(path).exists()


# ---------------------------------------------------------------------------
# Convenience accessors
# ---------------------------------------------------------------------------
def raw(name: str) -> Path:
    """Path to a large external input under data/raw/ (fetched by step00)."""
    return DATA_RAW / name


def cache(name: str) -> Path:
    """Path to a committed connector/LLM response under data/cache/."""
    return DATA_CACHE / name


def result(step: str, name: str) -> Path:
    """Path to a per-step output under results/<step>/."""
    d = STEP_DIRS.get(step, RESULTS / step)
    d.mkdir(parents=True, exist_ok=True)
    return d / name


def require(path: os.PathLike | str, hint: str = "") -> Path:
    """
    Assert that an input file exists, with an actionable error otherwise.

    Used at the top of each step to fail fast with a message pointing at the
    upstream step or fetch command that produces the missing file.
    """
    p = Path(path)
    if not p.exists():
        msg = f"Required input not found: {p}"
        if hint:
            msg += f"\n  -> {hint}"
        raise FileNotFoundError(msg)
    return p
