# data/raw/ — large external inputs (not committed)

This directory is populated by `src/step00_fetch_data.py`. Its contents (~8.6 GB
total) are git-ignored because they are large public datasets that should be
downloaded from source, not redistributed here:

| File | Source | Approx size |
|------|--------|-------------|
| `CRISPRGeneEffect.csv`, `Model.csv`, `CRISPRInferredCommonEssentials.csv`, `OmicsAbsoluteCNGene.csv`, `OmicsSomaticMutationsMatrix*.csv` | DepMap 24Q4 Public (figshare 27993248) | ~430 MB |
| `core_gbmap.h5ad` | GBmap core atlas (CELLxGENE Census) | ~8.1 GB |
| `GSE97930_*_snDrop-seq_UMI_Count_Matrix_*.txt.gz` | Lake et al. 2018 (GEO GSE97930) | ~47 MB |
| `TcgaTargetGTEX_phenotype.txt.gz` | UCSC Xena Toil | ~4 MB |

Run `python src/step00_fetch_data.py` to download them. In the default cached
mode the rest of the pipeline reproduces the published results from the committed
`data/cache/` and `results/` files **without** these raw inputs; they are only
required to recompute steps 1, 3, 4, 5, 6 from scratch.
