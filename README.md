# MuDataSeurat

## NOTE:

**This is a fork from PMBio's [MuDataSeurat](https://github.com/PMBio/MuDataSeurat)**

## This fork fixed issues with:

1. Compatible with Seurat v5
2. Tested and works with `anndata` (>=0.8), and `anndata-rs`(>=0.2.0).
3. Export all reductions, such as `UMAP`, `tSNE` etc. 
4. Fixed stack overflow issue because of obs column with NA
   - skip columns with all NA value
   - fixed string array with NA
4. Add two new keyword arguments to `WriteH5AD` and `WriteH5MU`: 
   - `scale.data`: wheter write `scale.data` to `anndata/mudata`
   - `sparse.type`: store `csc_matrix` or `csr_matrix` in `anndata/mudata`
 


Please ref to the main repo [MuDataSeurat](https://github.com/PMBio/MuDataSeurat) for more details
## Installation

Please install the main branch

```R
remotes::install_github("zqfang/MuDataSeurat")
```
