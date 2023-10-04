# MuDataSeurat

## NOTE:

**This is a fork from PMBio's [MuDataSeurat](https://github.com/PMBio/MuDataSeurat)**

## Features

1. Tested and works with `anndata` (>=0.8), and `anndata-rs`(>=0.2.0).
2. Export all reductions, such as `UMAP`, `tSNE` etc. 
3. Add two new keyword arguments to `WriteH5AD` and `WriteH5MU`: 
   - `scale.data`: wheter write `scale.data` to `anndata/mudata`
   - `sparse.type`: store `csc_matrix` or `csr_matrix` in `anndata/mudata`


Please ref to the main repo [MuDataSeurat](https://github.com/PMBio/MuDataSeurat) for more details
## Installation

```R
remotes::install_github("zqfang/MuDataSeurat")
```
