# Step 4 — Tumor-vs-TME Deconvolution (GBmap single-cell)

**Atlas:** Core GBmap (Ruiz-Moreno et al.), CELLxGENE collection 999f2a15-3d7e-440b-96ae-2c806799c08c,
dataset c888b684-6c51-431f-972a-6c963044cef0. **338,564 cells × 27,632 genes**, 17 annotated cell types.
**Expression:** log-normalized (X layer). **Compartments:** malignant cell (127,521; 37.7%) vs resolved TME
(macrophage, microglial, monocyte, dendritic, T, NK, B, oligodendrocyte, endothelial, mural, radial glial).

## Why this step

DepMap gives glioma-selective *dependency*, but a gene that is only expressed by tumor-associated
myeloid cells (GBM is 30–50% microglia/macrophage) is not a tumor-cell target even if the lineage screen
flags it. Single-cell separates malignant from microenvironment expression — the classic non-tumor-cell
expression trap (e.g. immune-checkpoint genes such as TIGIT/CTLA4 confined to the T-cell compartment, or
CD83 in the myeloid compartment — signals that come from distinct non-malignant compartments).

## Validation

Aggregation validated against canonical compartment markers: EGFR/SOX2/OLIG2 peak in malignant cells
(EGFR 0.94 malignant vs ~0 TME); CD68/AIF1 in macrophage/microglia; CD3E/PTPRC in T cells; PECAM1/VWF in
endothelial. Compartment means are correct. Malignant fraction 37.7% matches known GBM TME composition.

## Results

- **278 of 283 candidates measured** in GBmap (5 absent from the atlas gene set: KRTAP21-1, CSH2, TUBB,
  POTEI, NOPCHAP1 — recorded as not-measured, SC feature NaN, other features retained).
- **43 candidates flagged TME-dominated** (expression higher in a TME compartment than in malignant by
  >0.25 log units) — these are down-weighted: their dependency signal may reflect microenvironment biology,
  or they are poor tumor-cell targets.
- **PDGFRA is the top malignant-specific candidate** (specificity +0.31; malignant 0.36 vs max-TME 0.05) —
  the canonical GBM tumor-cell driver, cleanly tumor-restricted. Ideal profile.

### Most malignant-specific (top 10)
         sc_malignant_mean  sc_tme_max_mean      sc_tme_max_type  sc_malignant_specificity
PDGFRA               0.361            0.052           mural cell                     0.309
CCT6A                0.834            0.578           mural cell                     0.255
PDCD5                0.611            0.438           mural cell                     0.173
PPP1CB               0.703            0.533             monocyte                     0.169
POP7                 0.336            0.190      oligodendrocyte                     0.146
EIF3B                0.259            0.142      oligodendrocyte                     0.117
PPP2R1A              0.493            0.400           mural cell                     0.093
HSPA9                0.488            0.396      oligodendrocyte                     0.092
SETD5                0.322            0.243     endothelial cell                     0.079
NELFCD               0.260            0.186  natural killer cell                     0.074

### Most TME-dominated (bottom 10 — candidates to down-weight)
         sc_malignant_mean  sc_tme_max_mean    sc_tme_max_type  sc_malignant_specificity
RPS29                1.771            3.578             B cell                    -1.806
RPS17                1.269            3.010  radial glial cell                    -1.741
RPS3                 1.739            3.254             B cell                    -1.514
RPS25                1.821            3.175             B cell                    -1.354
RPL10A               1.568            2.838             B cell                    -1.271
RPS21                1.573            2.708  radial glial cell                    -1.135
RPS11                1.759            2.839             B cell                    -1.080
RPS9                 1.844            2.840  radial glial cell                    -0.996
ZFP36L1              0.615            1.609           monocyte                    -0.995
RPS16                2.061            3.026             B cell                    -0.965

## Honest interpretation & caveats

- **The negative signal is the reliable one.** `sc_tme_dominated_flag` robustly identifies candidates
  expressed in the wrong compartment (penalize). The positive `sc_malignant_specificity` score is weaker:
  in log-normalized space, ubiquitously/highly expressed housekeeping genes (ribosomal, chaperones —
  CCT6A, PPP1CB, C19orf53, ribosomal proteins) show a small positive malignant margin without being
  biologically "malignant-specific." Recommend `sc_malignant_specificity` as a **modest bonus**, not a
  strong driver; `sc_tme_dominated_flag` as a real down-weight.
- **Two specificity framings provided** (correlation 0.62): resolved-max-TME (`sc_malignant_specificity`)
  and collapsed neoplastic-vs-non-neoplastic (`sc_malignant_vs_tme_collapsed`). The collapsed metric is less
  sensitive to a single small TME population; both are in the table.
- **Ribosomal/housekeeping high-expressors** dominate the extreme-expression corner and are neither strong
  targets nor true TME traps — flagged in the figure for transparency.

## Features added to the matrix

`sc_malignant_mean`, `sc_tme_max_mean`, `sc_tme_max_type`, `sc_tme_mean`, `sc_pct_malignant_expressing`,
`sc_malignant_specificity`, `sc_tme_dominated_flag`, `sc_neoplastic_mean`, `sc_nonneoplastic_mean`,
`sc_malignant_vs_tme_collapsed`, `sc_measured`.

## Deliverables

- `singlecell_compartment_expression.csv` — per-candidate SC features (feeds the matrix)
- `gbmap_pseudobulk_mean.csv` / `gbmap_pseudobulk_pct.csv` — full 278×17 compartment matrices (checkpoint)
- `gbmap_celltype_counts.csv` — cells per compartment
- `singlecell_deconvolution.png` — malignant-vs-TME scatter + top-25 compartment heatmap
