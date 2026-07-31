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
5. Add two new keyword arguments to `WriteH5AD` and `WriteH5MU`: 
   - `scale.data`: whether write `scale.data` to `anndata/mudata` or not.
   - `sparse.type`: store `csc_matrix` or `csr_matrix` in `anndata/mudata`
6. Correctness fixes
   - `.h5mu` files are now readable by `mudata`: `obsm`, `varm`, `obsp`, `varp`,
     `obsmap` and `varmap` are always written.
   - Categorical (factor) columns keep their labels when a category is unused.
   - `NA` in string columns is stored as a missing value (as `categorical`,
     since anndata has no nullable-string encoding) instead of the text `"NaN"`.
   - `WriteH5AD` now errors on an unknown `assay` instead of silently writing
     a different one.
   - Reductions whose names contain non-alphanumeric characters (`umap_harmony`,
     `pca_uncorrected`, `mrVI_umap`, ...) are no longer dropped from `obsm`.
     Seurat strips those characters when it derives a reduction key
     (`umap_harmony` -> `umapharmony_`), which previously made such reductions
     fail the assay-matching check and be skipped silently. A reduction that is
     still not matched now raises a warning instead of disappearing.
7. Much faster on large objects. On a 500,000 cell object (2,000 features,
   150M nonzeros, 28 metadata columns):

   | | before | after |
   |---|---|---|
   | `ReadH5AD` | 54 s | 9 s |
   | `WriteH5AD` | 23 s | 23 s (`compression = "gzip"`, the default) |
   | `WriteH5AD` | 23 s | 3 s (`compression = "none"`) |

   - Reading a `csr_matrix` no longer converts it to CSC and then transposes it
     back: AnnData's CSR buffers already describe the matrix in the orientation
     Seurat wants, so it is built from them directly.
   - `nCount_*`/`nFeature_*` are no longer recomputed from the matrix when the
     stored `obs` already contains them (they were recomputed and then
     immediately overwritten). They are still computed when the file omits them.
   - `obs` and `var` are read once per file instead of once per consumer
     (`ReadH5AD` used to read `/obs` four times).
   - New `compression` argument on `WriteH5AD`/`WriteH5MU`: `"gzip"` (default,
     unchanged behaviour), `"none"`, or a gzip level from 0 to 9. Compression,
     not disk I/O, dominates write time; `"none"` is several times faster in
     exchange for a roughly 3x larger file.


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
seu <- ReadH5AD("export.h5ad")
seu <- ReadH5MU("export.h5mu")
```
You may also use the native support of anndata in R: `anndataR::read_h5ad`
