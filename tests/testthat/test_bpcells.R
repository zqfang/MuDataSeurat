library(Seurat)
library(MuDataSeurat)
library(Matrix)
library(hdf5r)
library(fs)  # for file_temp()

# Disk-backed (BPCells) matrices are written by streaming a block of cells at a
# time instead of materialising a dgCMatrix. The streaming writer is exercised
# here both directly -- with blocks served from an in-memory matrix, so the
# indptr arithmetic is covered without BPCells installed -- and end to end.

nobs <- 50
nvar <- 12

obs_names <- paste("obs", seq_len(nobs), sep = "-")
var_names <- paste("var", seq_len(nvar), sep = "-")

make_counts <- function() {
  set.seed(42)
  x <- rnbinom(n = nobs * nvar, prob = .9, size = 10)
  x <- Matrix(matrix(x, ncol = nobs), sparse = TRUE)  # => dgCMatrix, var x obs
  colnames(x) <- obs_names
  rownames(x) <- var_names
  x
}

# Read a csr_matrix group back the way the AnnData spec describes it, so the
# assertions do not depend on this package's own reader.
read_csr <- function(file, key = "X") {
  h5 <- H5File$new(file, mode = "r")
  on.exit(h5$close_all())
  grp <- h5[[key]]
  shape <- as.integer(h5attr(grp, "shape"))
  mat <- sparseMatrix(
    j = grp[["indices"]]$read(),
    p = grp[["indptr"]]$read(),
    x = grp[["data"]]$read(),
    dims = shape,
    index1 = FALSE
  )
  list(
    matrix = mat,
    shape = shape,
    encoding = h5attr(grp, "encoding-type")
  )
}

test_that("write_sparse_stream reproduces the matrix it was given in blocks", {
  counts <- make_counts()
  file <- paste0(file_temp(), ".h5")

  h5 <- H5File$new(file, mode = "w")
  # A block size that does not divide nobs evenly, so the ragged final block is
  # covered too.
  MuDataSeurat:::write_sparse_stream(
    h5, "X", dim(counts),
    function(j1, j2) counts[, j1:j2, drop = FALSE],
    block_cols = 7L
  )
  h5$close_all()

  out <- read_csr(file)
  expect_equal(out$encoding, "csr_matrix")
  expect_equal(out$shape, c(nobs, nvar))
  # The stream writes obs x var; the source is var x obs.
  expect_equal(unname(as.matrix(t(out$matrix))), unname(as.matrix(counts)))
})

test_that("write_sparse_stream handles a block size larger than the matrix", {
  counts <- make_counts()
  file <- paste0(file_temp(), ".h5")

  h5 <- H5File$new(file, mode = "w")
  MuDataSeurat:::write_sparse_stream(
    h5, "X", dim(counts),
    function(j1, j2) counts[, j1:j2, drop = FALSE],
    block_cols = 10000L
  )
  h5$close_all()

  expect_equal(unname(as.matrix(t(read_csr(file)$matrix))), unname(as.matrix(counts)))
})

test_that("write_sparse_stream writes an all-zero matrix as an empty csr_matrix", {
  counts <- make_counts()
  counts[] <- 0
  counts <- as(counts, "dgCMatrix")
  file <- paste0(file_temp(), ".h5")

  h5 <- H5File$new(file, mode = "w")
  MuDataSeurat:::write_sparse_stream(
    h5, "X", dim(counts),
    function(j1, j2) counts[, j1:j2, drop = FALSE],
    block_cols = 7L
  )
  h5$close_all()

  h5 <- H5File$new(file, mode = "r")
  on.exit(h5$close_all())
  # indptr still has one entry per observation plus the leading zero, and every
  # entry is zero because nothing was appended.
  expect_equal(length(h5[["X"]][["indptr"]]$read()), nobs + 1)
  expect_true(all(h5[["X"]][["indptr"]]$read() == 0))
  expect_equal(length(h5[["X"]][["data"]]$read()), 0)
})

test_that("indptr is written as a 64-bit integer", {
  counts <- make_counts()
  file <- paste0(file_temp(), ".h5")

  h5 <- H5File$new(file, mode = "w")
  MuDataSeurat:::write_sparse_stream(h5, "X", dim(counts),
                                     function(j1, j2) counts[, j1:j2, drop = FALSE])
  h5$close_all()

  h5 <- H5File$new(file, mode = "r")
  on.exit(h5$close_all())
  expect_equal(h5[["X"]][["indptr"]]$get_type()$to_text(), "H5T_STD_I64LE")
})

test_that("compression = \"none\" keeps the stream chunked but unfiltered", {
  counts <- make_counts()
  file <- paste0(file_temp(), ".h5")

  h5 <- H5File$new(file, mode = "w")
  MuDataSeurat:::write_sparse_stream(
    h5, "X", dim(counts),
    function(j1, j2) counts[, j1:j2, drop = FALSE],
    ds_args = MuDataSeurat:::resolve_compression("none")
  )
  h5$close_all()

  # An extendable dataset has to be chunked, so "none" can only turn the filter
  # off -- it cannot fall back to the contiguous layout the in-memory path uses.
  expect_equal(unname(as.matrix(t(read_csr(file)$matrix))), unname(as.matrix(counts)))
})

# ---------------------------------------------------------------------------
# End-to-end, with BPCells

skip_if_no_bpcells <- function() {
  skip_if_not_installed("BPCells")
}

make_bpcells_srt <- function(counts) {
  dir <- file_temp()
  mat <- BPCells::write_matrix_dir(
    BPCells::convert_matrix_type(counts, "uint32_t"),
    dir = dir
  )
  CreateSeuratObject(counts = mat)
}

test_that("a BPCells-backed assay round-trips through WriteH5AD", {
  skip_if_no_bpcells()

  counts <- make_counts()
  srt <- make_bpcells_srt(counts)
  expect_s4_class(SeuratObject::LayerData(srt, "counts"), "IterableMatrix")

  file <- paste0(file_temp(), ".h5ad")
  expect_true(WriteH5AD(srt, file))

  out <- read_csr(file)
  expect_equal(out$encoding, "csr_matrix")
  expect_equal(out$shape, c(nobs, nvar))
  expect_equal(unname(as.matrix(t(out$matrix))), unname(as.matrix(counts)))
})

test_that("a BPCells-backed assay round-trips back into a Seurat object", {
  skip_if_no_bpcells()

  counts <- make_counts()
  file <- paste0(file_temp(), ".h5ad")
  WriteH5AD(make_bpcells_srt(counts), file)

  srt <- ReadH5AD(file)
  expect_equal(dim(srt), c(nvar, nobs))
  expect_equal(rownames(srt), var_names)
  expect_equal(colnames(srt), obs_names)
  expect_equal(
    unname(as.matrix(SeuratObject::LayerData(srt, "counts"))),
    unname(as.matrix(counts))
  )
})

test_that("a lazily normalised BPCells matrix is written with the transform applied", {
  skip_if_no_bpcells()

  counts <- make_counts()
  srt <- NormalizeData(make_bpcells_srt(counts), verbose = FALSE)
  # NormalizeData leaves an unevaluated operation tree rather than a matrix.
  expect_s4_class(SeuratObject::LayerData(srt, "data"), "IterableMatrix")

  file <- paste0(file_temp(), ".h5ad")
  expect_true(WriteH5AD(srt, file))

  # counts and data are both present, so data becomes X (case 2).
  expect_equal(
    unname(as.matrix(t(read_csr(file, "X")$matrix))),
    unname(as.matrix(SeuratObject::LayerData(srt, "data")))
  )
  expect_equal(
    unname(as.matrix(t(read_csr(file, "layers/counts")$matrix))),
    unname(as.matrix(counts))
  )
})

test_that("a BPCells-backed assay round-trips through WriteH5MU", {
  skip_if_no_bpcells()

  counts <- make_counts()
  file <- paste0(file_temp(), ".h5mu")
  expect_true(WriteH5MU(make_bpcells_srt(counts), file))

  expect_equal(
    unname(as.matrix(t(read_csr(file, "mod/RNA/X")$matrix))),
    unname(as.matrix(counts))
  )
})

test_that("writing a disk-backed matrix as csc_matrix is refused", {
  skip_if_no_bpcells()

  file <- paste0(file_temp(), ".h5ad")
  expect_error(
    WriteH5AD(make_bpcells_srt(make_counts()), file, sparse.type = "csc_matrix"),
    "requires sparse.type"
  )
})

test_that("writing over the file a matrix is backed by is refused", {
  skip_if_no_bpcells()

  counts <- make_counts()
  file <- paste0(file_temp(), ".h5ad")
  WriteH5AD(make_bpcells_srt(counts), file)

  # Re-open the file we just wrote as the backing store, then try to write the
  # resulting object straight back over it.
  mat <- BPCells::open_matrix_anndata_hdf5(file)
  dimnames(mat) <- list(var_names, obs_names)
  srt <- CreateSeuratObject(counts = mat)

  expect_error(WriteH5AD(srt, file), "reads from that same location")
  expect_error(WriteH5MU(srt, file), "reads from that same location")
})

# ---------------------------------------------------------------------------
# Reading into a BPCells-backed object

# An .h5ad written from an ordinary in-memory object, to read back.
write_plain_h5ad <- function(counts, normalize = FALSE) {
  srt <- CreateSeuratObject(counts)
  if (normalize) {
    srt <- NormalizeData(srt, verbose = FALSE)
  }
  file <- paste0(file_temp(), ".h5ad")
  WriteH5AD(srt, file)
  file
}

test_that("the memory backend is the default and is unchanged", {
  counts <- make_counts()
  srt <- ReadH5AD(write_plain_h5ad(counts))

  expect_s4_class(LayerData(srt, "counts"), "dgCMatrix")
  expect_equal(unname(as.matrix(LayerData(srt, "counts"))), unname(as.matrix(counts)))
})

test_that("backend = 'bpcells' leaves the matrix on disk in the source file", {
  skip_if_no_bpcells()

  counts <- make_counts()
  file <- write_plain_h5ad(counts)
  srt <- ReadH5AD(file, backend = "bpcells")

  expect_true(inherits(LayerData(srt, "counts"), "IterableMatrix"))
  expect_equal(dim(srt), c(nvar, nobs))
  expect_equal(rownames(srt), var_names)
  expect_equal(colnames(srt), obs_names)
  expect_equal(unname(as.matrix(LayerData(srt, "counts"))), unname(as.matrix(counts)))
})

test_that("bpcells.dir converts each matrix into its own directory", {
  skip_if_no_bpcells()

  counts <- make_counts()
  file <- write_plain_h5ad(counts, normalize = TRUE)
  dir <- file_temp()
  srt <- ReadH5AD(file, backend = "bpcells", bpcells.dir = dir)

  expect_true(inherits(LayerData(srt, "counts"), "IterableMatrix"))
  # counts and data both present, so X is data and layers/counts is counts.
  expect_setequal(list.files(dir), c("X", "layers_counts"))
  expect_equal(unname(as.matrix(LayerData(srt, "counts"))), unname(as.matrix(counts)))
})

test_that("a converted object no longer depends on the source file", {
  skip_if_no_bpcells()

  counts <- make_counts()
  file <- write_plain_h5ad(counts)
  dir <- file_temp()
  srt <- ReadH5AD(file, backend = "bpcells", bpcells.dir = dir)

  # The whole point of converting: the object survives the .h5ad going away.
  unlink(file)
  expect_equal(unname(as.matrix(LayerData(srt, "counts"))), unname(as.matrix(counts)))
})

test_that("reductions and graphs stay in memory under the bpcells backend", {
  skip_if_no_bpcells()

  counts <- make_counts()
  srt <- CreateSeuratObject(counts)
  embeddings <- matrix(
    seq_len(nobs * 2) / 10, nrow = nobs,
    dimnames = list(obs_names, c("pca_1", "pca_2"))
  )
  srt[["pca"]] <- CreateDimReducObject(
    embeddings = embeddings, key = "pca_", assay = "RNA"
  )
  file <- paste0(file_temp(), ".h5ad")
  WriteH5AD(srt, file)

  back <- ReadH5AD(file, backend = "bpcells")
  # obsm is small and Seurat needs a real matrix for it.
  expect_true(is.matrix(Embeddings(back[["pca"]])))
  expect_equal(unname(Embeddings(back[["pca"]])), unname(embeddings))
})

test_that("ReadH5MU backs every modality and keeps their own feature names", {
  skip_if_no_bpcells()

  rna <- make_counts()
  adt <- make_counts()[seq_len(3), , drop = FALSE]
  rownames(adt) <- paste("adt", seq_len(3), sep = "-")

  srt <- CreateSeuratObject(rna, assay = "RNA")
  srt[["ADT"]] <- CreateAssay5Object(counts = adt)
  file <- paste0(file_temp(), ".h5mu")
  WriteH5MU(srt, file)

  dir <- file_temp()
  back <- ReadH5MU(file, backend = "bpcells", bpcells.dir = dir)

  expect_setequal(names(back@assays), c("RNA", "ADT"))
  for (assay in c("RNA", "ADT")) {
    expect_true(inherits(LayerData(back[[assay]], "counts"), "IterableMatrix"))
  }
  expect_equal(rownames(back[["ADT"]]), rownames(adt))
  expect_equal(
    unname(as.matrix(LayerData(back[["ADT"]], "counts"))),
    unname(as.matrix(adt))
  )

  # Each modality converts into its own subdirectory, and the names written
  # there are the modality's own -- BPCells would otherwise infer them from the
  # file's global /obs and /var, which for an .h5mu are the wrong names.
  expect_setequal(list.files(dir), c("RNA", "ADT"))
  on_disk <- BPCells::open_matrix_dir(file.path(dir, "ADT", "X"))
  expect_equal(rownames(on_disk), rownames(adt))
})

test_that("a disk-backed object read from a file can be written back out", {
  skip_if_no_bpcells()

  counts <- make_counts()
  file <- write_plain_h5ad(counts)
  srt <- ReadH5AD(file, backend = "bpcells")

  out <- paste0(file_temp(), ".h5ad")
  expect_true(WriteH5AD(srt, out))
  expect_equal(unname(as.matrix(t(read_csr(out)$matrix))), unname(as.matrix(counts)))
})

test_that("backend = 'bpcells' rejects bpcells.dir on the memory backend", {
  file <- write_plain_h5ad(make_counts())
  expect_error(ReadH5AD(file, bpcells.dir = "somewhere"), "only meaningful")
  expect_error(
    ReadH5MU(paste0(file_temp(), ".h5mu"), bpcells.dir = "somewhere"),
    "only meaningful"
  )
})

test_that("an unknown backend is rejected", {
  file <- write_plain_h5ad(make_counts())
  expect_error(ReadH5AD(file, backend = "hdf5array"))
})
