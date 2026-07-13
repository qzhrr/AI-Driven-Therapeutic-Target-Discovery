"""
kernel.py — reusable scoring primitives for the cancer target-discovery pipeline.

Auto-loaded into the live kernel when this skill is loaded. These are the
cancer-independent building blocks of the 12-step method (see SKILL.md), factored
out so they can be reused when adapting the pipeline to a new cancer type or lane.

All functions are pure and operate on pandas Series / numpy arrays. Heavy imports
are deferred into function bodies so loading the skill is cheap.

Available after load:
    cohens_d(group_a, group_b)                  effect size (pooled SD)
    benjamini_hochberg(pvalues)                 BH-FDR adjusted p-values
    welch_ttest(group_a, group_b)               Welch t statistic + p per row
    minmax01(series, zero_preserving=False)     scale to [0, 1]
    graded_percentile_safety(series, log=True)  1 - percentile rank (higher = safer)
    weighted_composite(category_scores, weights)  weighted sum across categories
    classify_masked_effect(freq_named, freq_masked)  evidence vs name-recognition class
    masked_ablation_effect(df, named, masked)   add an 'effect' column to a table
"""

# Category order used throughout the pipeline (kept as a literal constant).
EVIDENCE_CATEGORIES = [
    "dependency",
    "selectivity",
    "synthetic_lethality",
    "tumor_specificity",
    "tumor_enrichment",
    "therapeutic_window",
    "tractability",
]


def cohens_d(group_a, group_b):
    """
    Cohen's d effect size between two 1-D samples, using the pooled standard
    deviation with an (n-2) denominator (as used in the dependency step).

    Positive d means group_a > group_b. For DepMap gene-effect values (lower =
    more essential), pass glioma as group_a and non-glioma as group_b, then use
    -d as the selectivity signal.
    """
    import numpy as np

    a = np.asarray(group_a, dtype=float)
    b = np.asarray(group_b, dtype=float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    n_a, n_b = len(a), len(b)
    if n_a < 2 or n_b < 2:
        return float("nan")
    pooled_sd = np.sqrt(((n_a - 1) * a.std(ddof=1) ** 2 + (n_b - 1) * b.std(ddof=1) ** 2) / (n_a + n_b - 2))
    if pooled_sd == 0:
        return 0.0
    return float((a.mean() - b.mean()) / pooled_sd)


def welch_ttest(group_a, group_b):
    """Welch's t-test (unequal variance). Returns (t_statistic, p_value)."""
    from scipy import stats

    t, p = stats.ttest_ind(group_a, group_b, equal_var=False, nan_policy="omit")
    return float(t), float(p)


def benjamini_hochberg(pvalues):
    """Benjamini-Hochberg FDR-adjusted p-values (returns a numpy array)."""
    from statsmodels.stats.multitest import multipletests
    import numpy as np

    p = np.asarray(pvalues, dtype=float)
    mask = ~np.isnan(p)
    out = np.full_like(p, np.nan, dtype=float)
    if mask.sum() > 0:
        out[mask] = multipletests(p[mask], method="fdr_bh")[1]
    return out


def minmax01(series, zero_preserving=False):
    """
    Scale a series to [0, 1]. If zero_preserving, a series with no spread maps to
    all zeros (used for synthetic-lethality features where 0 means 'no SL'),
    matching the pipeline's mm01 helper.
    """
    import pandas as pd

    s = pd.Series(series).astype(float)
    lo, hi = s.min(), s.max()
    if hi > lo:
        return (s - lo) / (hi - lo)
    return s * 0 if zero_preserving else s.clip(lower=0, upper=0)


def graded_percentile_safety(series, log=True):
    """
    Graded safety score: 1 - percentile-rank of (optionally log1p) expression.
    Higher = expressed lower across the reference tissues = safer / wider window.
    This is the pipeline's therapeutic-window scoring (step 5).
    """
    import numpy as np
    import pandas as pd

    s = pd.Series(series).astype(float).copy()
    if log:
        s = np.log1p(s)
    return 1.0 - s.rank(pct=True)


def weighted_composite(category_scores, weights):
    """
    Weighted sum of per-category scores.

    category_scores: DataFrame (index = gene, columns = category names in [0,1]).
    weights: dict or Series mapping category -> weight (need not be normalized;
             they are normalized to sum 1 here).
    Returns a Series of composite scores indexed by gene.
    """
    import pandas as pd

    scores = pd.DataFrame(category_scores)
    w = pd.Series(weights, dtype=float)
    cats = [c for c in scores.columns if c in w.index]
    if not cats:
        raise ValueError("no overlapping categories between scores and weights")
    wn = w[cats] / w[cats].sum()
    return (scores[cats] * wn).sum(axis=1)


def classify_masked_effect(freq_named, freq_masked):
    """
    Classify one candidate by how its LLM nomination frequency changes when the
    gene name is masked (the step-11a bias test):

      'evidence-carried (robust)'  named>=0.5 and masked>=0.5*named
      'NAME-INFLATED (drops when masked)'  named>=0.5 and masked<0.5*named
      'NAME-SUPPRESSED (rises when masked)'  named<0.3 and masked>=0.5
      'other/low' otherwise
    """
    if freq_named >= 0.5 and freq_masked >= 0.5 * freq_named:
        return "evidence-carried (robust)"
    if freq_named >= 0.5 and freq_masked < 0.5 * freq_named:
        return "NAME-INFLATED (drops when masked)"
    if freq_named < 0.3 and freq_masked >= 0.5:
        return "NAME-SUPPRESSED (rises when masked)"
    return "other/low"


def masked_ablation_effect(df, named_col="freq_named", masked_col="freq_masked",
                           active_threshold=0.10):
    """
    Add 'delta' and 'effect' columns to a named-vs-masked nomination table.

    df must have the two frequency columns. Genes below active_threshold in both
    arms are labeled 'not-nominated'. Returns a copy sorted by named frequency.
    """
    import pandas as pd

    out = pd.DataFrame(df).copy()
    out["delta"] = out[masked_col] - out[named_col]
    active = (out[named_col] >= active_threshold) | (out[masked_col] >= active_threshold)
    out["effect"] = "not-nominated"
    out.loc[active, "effect"] = [
        classify_masked_effect(n, m)
        for n, m in zip(out.loc[active, named_col], out.loc[active, masked_col])
    ]
    return out.sort_values(named_col, ascending=False)
