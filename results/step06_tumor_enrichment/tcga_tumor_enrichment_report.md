# Step 6 — Tumor expression enrichment (TCGA-GBM vs GTEx normal brain)

## Purpose
Add the bulk patient-cohort over-expression axis: is each candidate up-regulated in GBM tumor relative to
normal brain, and in what fraction of patients? Complements the Step-4 single-cell signal ("expressed in tumor
*cells*") with cohort-level "over-expressed across *patients*".

## Data source — UCSC Xena Toil hub (uniform pipeline, no batch confound)
The critical design choice: TCGA-GBM tumor and GTEx normal brain are both taken from the **Xena Toil
`TcgaTargetGtex_rsem_gene_tpm`** dataset — TCGA and GTEx reprocessed through **one identical RNA-seq pipeline**
(values = log2(TPM+0.001)). This removes the cross-pipeline batch confound that would contaminate a raw
TCGA-GDC-vs-GTEx comparison. Retrieved by streaming the 1.3 GB matrix once and keeping only the 283 candidate
gene rows (285 Ensembl IDs; 11 recent HUGO renames recovered via GENCODE-v23 legacy symbols → 283/283 covered).

- **TCGA-GBM primary tumor:** 153 samples
- **GTEx normal brain:** 1,152 samples (13 subregions)

## Features (per candidate)
- `tcga_gbm_median_log2tpm`, `gtex_brain_median_log2tpm`
- `tcga_vs_normal_log2fc` = median tumor − median normal (log2 fold-change) — core over-expression signal
- `tcga_pct_over_normal` = fraction of GBM tumors above the normal-brain 75th percentile (prevalence)
- `tcga_vs_normal_p` / `tcga_vs_normal_fdr` (Mann-Whitney U + BH)

## Results
- **Distribution is discriminating, not inflated:** 169 over-expressed (log2FC>1), 110 neutral, 4 under-expressed
  (median +1.14). Under-expressed tail (AMY2A, GPR61, HJV) makes biological sense — non-brain genes.
- **Top tumor-enriched are the GBM proliferation signature:** UBE2C, SGO1, KIF18B, SKA1, KIF18A, KIF2C, ERCC6L, H3C8 — cell-cycle/mitotic
  genes, expected to be up vs post-mitotic normal brain.
- **Key targets validate:** MET (+2.16, 63% of tumors), PDGFRA (+1.74, 66%), WDR77 (+1.76, 99%), TEAD1 (+1.75, 88%),
  E2F3 (+1.57, 81%) — all clearly tumor-enriched.
- **PTK2/FAK is neutral (+0.13):** a dependency that is NOT transcriptionally over-expressed — consistent with its
  high *normal*-brain expression (Step 5). Good internal consistency.

## Bulk vs single-cell: orthogonal, by design (r=0.01)
Panel B shows bulk tumor-enrichment (log2FC) and single-cell malignant expression are essentially uncorrelated.
This is expected and useful, NOT a discrepancy:
- **log2FC** is a *contrast* vs normal brain — high for anything a post-mitotic brain lacks (cell-cycle genes light up).
- **sc_malignant_mean** is an *absolute level* in malignant cells — high for housekeeping genes regardless of tumor-specificity.
They capture different information, so they are **non-redundant features** for the matrix. A candidate strong on
both (e.g. WDR77, MET, PDGFRA) has convergent bulk + single-cell support; one strong on only one axis is flagged
by their disagreement.

## Deliverables
- `tcga_tumor_enrichment.csv` — 283 candidates × 6 features (feeds the matrix)
- `tcga_tumor_enrichment.png` — tumor-vs-normal volcano + bulk-vs-single-cell cross-check
