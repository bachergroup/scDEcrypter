# scDEcrypter: Uncertainty-aware differential expression analysis for viral infection in scRNA-seq
scDEcrypter models infection status and other cell-level variables as latent variables with partial observability. scDEcrypter implements a regularized two-way mixture model, where mixture weights estimate cells' probabilistic membership to combinations of cell states (e.g., infection status and cell type). The resulting weights are used to estimate cell-state-specific mean expression profiles and to account for cell-state uncertainty in differential expression testing.

Additional details and a FAQ for scDEcrypter are described in the vignette: https://github.com/bachergroup/scDEcrypter/blob/main/vignettes/scDEcrypter.pdf

## Installation
```R
library(devtools)
devtools::install_github("https://github.com/bachergroup/scDEcrypter")
library(scDEcrypter)
```

## Author
Luer Zhong <luerzhong@ufl.edu>

Rhonda Bacher <rbacher@ufl.edu>

## Cite

