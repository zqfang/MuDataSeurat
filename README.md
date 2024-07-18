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

## Usage


### Export to H5AD
```R
# MuDataSeurat only export 3 layers: count, data, scale.data
# so, for seurat v5, join layer first, for each modality

library(MuDataSeurat)

DefaultAssay(seu) = "RNA"
seu = JoinLayers(seu)
WriteH5AD(seu, "export.h5ad",  assay="RNA", scale.data = FALSE, overwrite=TRUE)

## h5mu file for multi-modality

WriteH5MU(seu, "export.h5mu", overwrite=TRUE)
```

### Read H5AD to Seurat

```R
ReadH5AD()
ReadH5MU()
```
You may use native support of anndata, via `anndataR::read_h5ad`
