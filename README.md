# scDEcrypter: Uncertainty-aware differential expression analysis for viral infection in scRNA-seq

scDEcrypter is an R package for differential expression (DE) analysis of single-cell RNA-seq data in which infection status, and optionally other cell-level labels such as cell type, are only partially known. It fits a penalized two-way mixture model in which a small set of confidently labeled cells anchor the latent states. The model estimates each cell's probability of belonging to each combination of states (e.g., infection status and cell type).

The cells are split into a generation set and a test set. The generation set is used to fit the model. In the test set, the inferred state probabilities (weights) are used to estimate state-specific mean expression and to test for infection-associated genes within each cell type while accounting for state uncertainty with p-values from a permutation test.

A step-by-step tutorial, including key assumptions and guidance on parameter and gene selection, is in the vignette: https://github.com/bachergroup/scDEcrypter/tree/main/vignettes

## Installation

The Bioconductor package `transformGamPoi` is required and is not installed automatically.

```R
install.packages(c("BiocManager", "devtools"))
BiocManager::install("transformGamPoi")
devtools::install_github("bachergroup/scDEcrypter")
library(scDEcrypter)
```

## Authors

- Luer Zhong <luerzhong@ufl.edu>
- Aaron Molstad <amolstad@umn.edu>
- Rhonda Bacher <rbacher@ufl.edu>

## Citation

If you use scDEcrypter, please cite the bioRxiv pre-print

Luer Zhong, Karl Ensberg, Scott Tibbets, Aaron J. Molstad, Rhonda Bacher. scDEcrypter: Uncertainty-aware differential expression analysis for viral infection in scRNA-seq. 
bioRxiv (pre-print). doi: https://doi.org/10.64898/2026.03.09.710583 

The release used in the manuscript is archived at Zenodo: https://doi.org/10.5281/zenodo.23093820