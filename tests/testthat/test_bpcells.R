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
