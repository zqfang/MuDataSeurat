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
5. Add one new keyword arguments to `WriteH5AD` and `WriteH5MU`:
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
   - `scale.data` is no longer exported, and the `scale.data` argument is
     deprecated (still accepted, ignored, and warned about only when `TRUE`).
     `ScaleData()` runs on the variable features only while AnnData requires `X`
     to span every feature, so exporting it meant padding the dense scaled
     matrix out with `NaN` to `n_features / n_variable_features` times its size
     — commonly 10-15x, and the largest allocation a write ever made. The
     resulting `X` was mostly `NaN`, which is not something a reader can compute
     on. Recompute scaling after reading instead.
7. Much faster on large objects. On a 500,000 cell object (2,000 features,
   150M nonzeros, 28 metadata columns):

   |               | before | after                                        |
   | ------------- | ------ | -------------------------------------------- |
   | `ReadH5AD`  | 54 s   | 9 s                                          |
   | `WriteH5AD` | 23 s   | 23 s (`compression = "gzip"`, the default) |
   | `WriteH5AD` | 23 s   | 3 s (`compression = "none"`)               |


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
8. **Disk-backed ([BPCells](https://github.com/bnprks/BPCells)) matrices, both
   directions.** `ReadH5AD`/`ReadH5MU` can leave the count matrices on disk
   instead of loading them (`backend = "bpcells"`), and `WriteH5AD`/`WriteH5MU`
   stream such an object back out a block of cells at a time. Neither direction
   ever holds the whole matrix in memory. `ConvertToSeuratBPCells()` and
   `ConvertToSeuratInMemory()` move an object you already have between the
   two. See
   [Working with BPCells (on-disk matrices)](#working-with-bpcells-on-disk-matrices).
9. Reading now always produces a **Seurat v5 assay** (`Assay5`), and several
   metadata bugs are fixed along with it:

   - Feature metadata on a v5 assay used to be dropped entirely on write. Since
     `CreateSeuratObject` has returned v5 assays by default since Seurat 5, this
     meant `/var` was written empty for almost everyone.
   - `/var` columns were duplicated on read (attached once when the assay was
     built, then `cbind`ed on again).
   - Writing an integer column containing `NA`s recursed until the stack
     overflowed. `FindVariableFeatures` produces such a column
     (`vf_vst_counts_rank`), so this was easy to hit.
   - Logical columns are stored as a real two-value boolean. hdf5r's default
     three-value enum (`FALSE`/`TRUE`/`NA`) is read as `uint8` by `h5py`, which
     made `highly_variable` come out as an integer and made columns with `NA`s
     unreadable by `anndata` outright.

## Installation

Please install the main branch

```R
remotes::install_github("zqfang/MuDataSeurat")
```

## Usage

### Export to H5AD, H5MU

MuDataSeurat exports 2 layers: `counts` and `data`. When both are present,
`counts` becomes `layers['counts']` and `data` becomes `X`; when only one is
present it becomes `X`.

`scale.data` is **not** exported. `ScaleData()` runs on the variable features
only, while AnnData requires `X` to span every feature, so writing it meant
padding the dense scaled matrix out to the full var axis with `NaN` — many times
the size of the real data, and mostly `NaN` once written. Recompute it after
reading instead (`ScaleData()` in Seurat, `sc.pp.scale` in scanpy).

You need `JoinLayers` for each modality first with seurat v5.

```R
library(MuDataSeurat)

DefaultAssay(seu) = "RNA"
seu = JoinLayers(seu) # critical for seurat v5

## write unimodal h5ad
WriteH5AD(seu, "export.h5ad",  assay="RNA", overwrite=TRUE)

## write multimodal h5mu
WriteH5MU(seu, "export.h5mu", overwrite=TRUE)
```

### Read H5AD to Seurat

```R
seu <- ReadH5AD("export.h5ad")
seu <- ReadH5MU("export.h5mu")
```

Both return a Seurat v5 object: matrices become v5 layers (`counts`, `data`,
`scale.data`), and `/var` becomes the assay's feature metadata, reachable with
`seu[["RNA"]][[]]`. For files too large to fit in memory, see
[Working with BPCells](#working-with-bpcells-on-disk-matrices).

### Working with BPCells (on-disk matrices)

[BPCells](https://github.com/bnprks/BPCells) keeps count matrices on disk and
streams them, so an object larger than memory can still be worked with. Both
directions are supported: reading can leave the matrices on disk, and writing can
stream them back out.

BPCells is an optional dependency; install it separately:

```R
remotes::install_github("bnprks/BPCells/r")
```

#### Reading into a BPCells-backed object

Pass `backend = "bpcells"` to keep `X`, `layers` and `raw` on disk:

```R
# matrices stay in the .h5ad; nothing is copied
seu <- ReadH5AD("big.h5ad", backend = "bpcells")

# or convert them once into BPCells' own format
seu <- ReadH5AD("big.h5ad", backend = "bpcells", bpcells.dir = "big_bpcells")

seu <- ReadH5MU("big.h5mu", backend = "bpcells", bpcells.dir = "big_bpcells")
```

The two forms trade off differently:

| `bpcells.dir` | cost to open | afterwards |
| ------------- | ------------ | ---------- |
| `NULL` (default) | nothing is copied | every pass re-reads the HDF5; the object breaks if the source file moves or is deleted |
| a directory | one pass over the data, plus disk | much faster to compute on, and the object survives `saveRDS` and the source file going away |

Reading a 120,000 cell file (2,000 features, 391 MB uncompressed `.h5ad`):

| backend | peak memory | time  |
| ------- | ----------- | ----- |
| `"memory"` (default) | 858 MB | 2.3 s |
| `"bpcells"`, no dir  | 145 MB | 2.2 s |
| `"bpcells"` + dir    | 149 MB | 1.9 s |

Only the count matrices are affected. `obsm`, `varm` and `obsp` — reductions and
graphs — are always read into memory, since Seurat needs real matrices for them
and they are far smaller.

With `.h5mu`, each modality converts into its own subdirectory
(`big_bpcells/RNA/X`, `big_bpcells/ADT/X`, ...) carrying that modality's own
feature names.

#### Writing a BPCells-backed object

`WriteH5AD`/`WriteH5MU` accept such objects directly and stream them out a block
of cells at a time, rather than pulling the whole matrix into memory first. This
works whether the object came from `ReadH5AD(..., backend = "bpcells")` or was
built with BPCells directly:

```R
library(BPCells)
library(MuDataSeurat)

# One-time conversion of the counts to BPCells' on-disk format
mat <- open_matrix_anndata_hdf5("big.h5ad")   # or open_matrix_10x_hdf5(), etc.
mat <- convert_matrix_type(mat, "uint32_t")   # raw counts only -- see below
mat <- write_matrix_dir(mat, dir = "big_bpcells")

seu <- CreateSeuratObject(counts = open_matrix_dir("big_bpcells"))
seu <- NormalizeData(seu)

WriteH5AD(seu, "export.h5ad")
WriteH5MU(seu, "export.h5mu")
```

`convert_matrix_type(mat, "uint32_t")` matters more than it looks: BPCells'
bitpacking compresses integers far better than floats, and skipping it can leave
the on-disk directory several times larger than it needs to be (5.4x on a small
test matrix). Use it only for raw counts — converting a normalized matrix would
truncate it to integers. Note that this package always writes matrix values as
float64, so counts read back out of an `.h5ad` arrive as doubles even when they
are whole numbers; BPCells will warn about poor compression if you skip the
conversion.

Normalized layers work too. `NormalizeData` on a disk-backed assay leaves an
unevaluated operation tree rather than a matrix, and the transform is applied as
each block is pulled, so nothing is materialised in full.

On a 120,000 cell object (2,000 features, 63.6M nonzeros), writing the same data
from a BPCells-backed assay versus an in-memory one:

|                     | peak memory | time  |
| ------------------- | ----------- | ----- |
| BPCells (streamed)  | 213 MB      | 3.0 s |
| in-memory `dgCMatrix` | 946 MB    | 1.7 s |

The output files are byte-for-byte identical apart from the `indptr` integer
width. Streaming trades some speed for a bounded memory ceiling.

Two constraints to be aware of:

- **`sparse.type` must be `"csr_matrix"`** (the default). A `csc_matrix` is
  grouped by feature, which is the opposite of the order a disk-backed matrix
  streams in; producing one would mean rewriting the whole matrix with
  `BPCells::transpose_storage_order()` first, so it is refused rather than done
  silently.
- **You cannot write over the file a matrix reads from.** A BPCells matrix keeps
  reading from its source for as long as the object is alive, so
  `WriteH5AD(seu, "big.h5ad")` where `seu`'s counts are backed by `big.h5ad` is
  refused — it would destroy the data being read. Write elsewhere.

The two can be chained, so a large file can be re-encoded without ever holding
its matrix in memory:

```R
seu <- ReadH5AD("big.h5ad", backend = "bpcells")
WriteH5MU(seu, "big.h5mu")   # streamed straight back out
```

#### Converting an object you already have

Two helpers move an existing Seurat object between memory and disk, leaving
everything else — feature metadata, variable features, reductions, graphs —
untouched:

```R
# in-memory -> disk-backed
seu <- ConvertToSeuratBPCells(seu, dir = "seu_bpcells")

# disk-backed -> in-memory
seu <- ConvertToSeuratInMemory(seu)
```

`ConvertToSeuratBPCells` writes each layer to `<dir>/<assay>/<layer>`, so the
object becomes self-contained and survives `saveRDS`. It also works on an object
read with `backend = "bpcells"` and no `bpcells.dir`, which is how you detach
such an object from the `.h5ad` it is still reading from.

Two defaults worth knowing:

- **`scale.data` stays in memory.** It is dense, usually covers only the variable
  features, and gains little from being moved. Pass `layers = "scale.data"` to
  convert it anyway.
- **No type casting happens unless you ask.** Casting raw counts to `"uint32_t"`
  compresses much better, but the same cast would truncate a normalized layer, so
  name the layers it applies to:

  ```R
  seu <- ConvertToSeuratBPCells(seu, dir = "seu_bpcells",
                                type = c(counts = "uint32_t"))
  ```

  Passing a bare `type = "uint32_t"` applies to every layer; if any in-memory
  layer does not hold non-negative integers, the call is refused rather than
  silently corrupting it.

`ConvertToSeuratInMemory` materializes the matrices, so check `dim(seu)` first —
that is the point of the call, but it does mean the object has to fit.
