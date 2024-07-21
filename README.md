# MuDataSeurat

**This is a fork from PMBio's [MuDataSeurat](https://github.com/PMBio/MuDataSeurat)**

Please refer to the original repo for more details

## Why this fork ?

I find `MuDataSeurat` to be the most compatible tool for converting `SeuratObject` to `H5AD/H5MU` format. It does not require python runtime as well.  

The original repository activity seems quite low, and unfortunately, the bugs have not been fixed promptly. Since I use `MuDataSeruat` quite often, I've decided to create my own fork and make my own version. I will do my best to ensure that it remains compatible with the latest pull requests from the original repository.

## New features

1. Compatible with Seurat v5
2. Tested and works with `anndata` (>=0.8), and `anndata-rs`(>=0.2.0).
3. Export all missing reductions, such as `UMAP`, `tSNE` etc. 
4. Fixed stack overflow issue because of obs column containing NAs
   - skip columns with all NA value
   - fixed string array with NA
4. Add two new keyword arguments to `WriteH5AD` and `WriteH5MU`: 
   - `scale.data`: whether write `scale.data` to `anndata/mudata` or not.
   - `sparse.type`: store `csc_matrix` or `csr_matrix` in `anndata/mudata`


## Installation

Please install the main branch

```R
remotes::install_github("zqfang/MuDataSeurat")
```

## Usage


### Export to H5AD, H5MU

MuDataSeurat only export 3 layers: `count`, `data`, `scale.data`

Therefore, You need `JoinLayers` for each modality first with seurat v5.

```R
library(MuDataSeurat)

DefaultAssay(seu) = "RNA"
seu = JoinLayers(seu) # critical for seurat v5

## write unimodal h5ad
WriteH5AD(seu, "export.h5ad",  assay="RNA", scale.data = FALSE, overwrite=TRUE)

## write multimodal h5mu
WriteH5MU(seu, "export.h5mu", overwrite=TRUE)
```

### Read H5AD to Seurat

```R
ReadH5AD()
ReadH5MU()
```
You may also use the native support of anndata in R: `anndataR::read_h5ad`
