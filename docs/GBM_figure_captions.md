# Nature-style figure captions — GBM target-discovery pipeline (v4)

All statistics recomputed from raw source tables (DepMap 24Q4; final_candidate_ranking.csv, 283 candidates;
gene_masked_ablation.csv; kif2c_glioma_dependency_perline.csv; glioma_dependency_scored.parquet).
Per-target accent colours are consistent across all figures: KIF2C = vermillion, VRK1 = blue, ELAVL1 = green;
KIF18A = grey (counter-example); positive controls = orange.

---

**Figure 1 | The pipeline recovers known GBM biology from a whole-genome start.**
(a) Prioritisation funnel: a whole-genome CRISPR dependency screen (17,916 genes, DepMap 24Q4) is reduced to
283 glioma-selective dependency candidates, scored on a 7-category composite augmented by an expression-context
synthetic-lethality layer (SL-v2), then triaged to three leads (KIF2C, VRK1, ELAVL1) by a blinded gene-masking
ablation and prior-art review. (b) Recovery of declared dependency positive controls among the 283 candidates
(composite rank, log scale). Four of five controls fall in the top 12 (PDGFRA #2, MET #4, VRK1 #7, FGFR1 #12);
PTPN11 (#74) is a pan-essential gene and is correctly down-weighted by the selectivity-aware composite. Controls
rank higher than non-controls overall (one-sided Mann–Whitney U, p = 3.1×10⁻⁵). Source: final_candidate_ranking.csv
(rank_all, is_dependency_poscontrol).

**Figure 2 | Dependency-first screening surfaces a 9p21-deletion synthetic-lethal lead (KIF2C).**
(a) Glioma-selective dependency landscape: Cohen's d (glioma vs all other lineages) versus −log₁₀ FDR for all
17,916 genes; the 283 candidates are highlighted. Portfolio leads and positive controls are labelled. EGFR sits
at positive d (rank 17,834) — glioma lines are *less* dependent on EGFR after copy-number correction — and is
excluded as a CRISPR/amplification artefact. Source: glioma_dependency_scored.parquet. (b) KIF2C dependency
stratified by 9p21/CDKN2A status across glioma lines: CDKN2A-deleted lines are more dependent on KIF2C than
CDKN2A-intact lines (mean gene effect −0.503 vs −0.353; Cohen's d = −0.70; two-sided Mann–Whitney p = 8.3×10⁻³;
n = 44 deleted / 28 intact). KIF18A shows no genotype selectivity (d = −0.09, n.s.; uniformly essential),
distinguishing a synthetic-lethal dependency from a pan-essential one. Source: kif2c_glioma_dependency_perline.csv.

**Figure 3 | Multi-evidence feature matrix for the top-25 candidates.**
Category scores (0–1) across the seven evidence categories (dependency, selectivity, synthetic lethality,
tumour specificity, tumour enrichment, therapeutic window, tractability), rows ordered by final composite rank
(rank_all). Portfolio leads are coloured (KIF2C, VRK1, ELAVL1); KIF18A is grey; dependency positive controls are
olive. The right-hand bar encodes novelty tier (approved small-molecule drug / has small-molecule chemistry /
novel-undrugged). Source: feature_matrix_scored.csv.

**Figure 4 | A blinded gene-masking test separates evidence from name recognition.**
(a) Nomination frequency (100 independent runs) for candidates presented with their real gene symbol (open marker)
versus anonymised opaque codes (filled marker); scores were byte-identical between arms, so only the name changed.
Genes are coloured by effect class: evidence-carried (blue; frequency is preserved when masked — VRK1 1.00→0.99,
ELAVL1 0.98→0.92, MET, PDGFRA), name-inflated (orange; frequency collapses when masked — KIF18A 1.00→0.00,
CDK7 1.00→0.03, KIF2C 0.86→0.29), and name-suppressed (pink; JUN 0.03→1.00, nominated only once its oncogene name
is hidden). KIF18A's collapse to zero is why it was removed from the headline. (b) Masked nomination frequency
tracks the pipeline's own selectivity evidence far better than named frequency (Spearman ρ = 0.12 named → 0.57
masked), i.e. removing names re-aligns the model with the data. 87 of 100 masked runs returned parseable output.
Source: gene_masked_ablation.csv.

**Figure 5 | The triage axes: therapeutic window and druggability separate the three leads.**
(a) Safety composite (higher = safer normal-tissue profile) versus tumour enrichment (TCGA-GBM vs normal,
log₂ fold change) across the 283 candidates. KIF2C sits in the safe-and-enriched quadrant; VRK1 is enriched with a
moderate safety profile; ELAVL1 is enriched but lower-safety. (b) All three leads are undrugged (novelty score = 1.0),
so they are separated by classical small-molecule tractability: VRK1 (0.80) is the cleanest small-molecule pick,
KIF2C (0.60) matches the KIF18A reference, and ELAVL1 (0.55) is lowest for classical inhibition — but ELAVL1/HuR has
an active targeted-degrader chemistry route that the classical tractability score does not capture. Source:
final_candidate_ranking.csv (safety_composite, tcga_vs_normal_log2fc, tractability_score, novelty_score).
