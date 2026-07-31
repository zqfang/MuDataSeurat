library(Seurat)
library(MuDataSeurat)
library(Matrix)
library(hdf5r)
library(fs) # for file_temp()

# Sparse matrix storage round-trips.
# Seurat stores features x cells; AnnData stores obs x var, so the shape
# attribute is reversed on write. The default storage type is csr_matrix.

nobs <- 10
nvar <- 20

obs_names <- paste("obs", 1:nobs, sep = "-")
var_names <- paste("var", 1:nvar, sep = "-")

# Hard-coded value to inject
true_val <- 0.1234569
true_val_i <- 3
true_val_j <- 7

make_srt <- function() {
  x <- rnbinom(n = nobs * nvar, prob = .95, size = 10)
  x <- Matrix(matrix(x, ncol = nobs), sparse = TRUE) # => dgCMatrix
  x[true_val_i, true_val_j] <- true_val
  colnames(x) <- obs_names
  rownames(x) <- var_names
  CreateSeuratObject(counts = x)
}

fileh5ad_r <- paste0(file_temp(), ".h5ad")
fileh5ad_c <- paste0(file_temp(), ".h5ad")

test_that("dgCMatrix is written to .h5ad as csr_matrix by default", {
  expect_true(WriteH5AD(make_srt(), fileh5ad_r))

  h5 <- H5File$new(fileh5ad_r, mode = "r")
  on.exit(h5$close_all())
  expect_equal(h5attr(h5[["X"]], "encoding-type"), "csr_matrix")
  expect_equal(as.integer(h5attr(h5[["X"]], "shape")), c(nobs, nvar))
})

test_that("dgCMatrix can be written to .h5ad as csc_matrix", {
  expect_true(WriteH5AD(make_srt(), fileh5ad_c, sparse.type = "csc_matrix"))

  h5 <- H5File$new(fileh5ad_c, mode = "r")
  on.exit(h5$close_all())
  expect_equal(h5attr(h5[["X"]], "encoding-type"), "csc_matrix")
  expect_equal(as.integer(h5attr(h5[["X"]], "shape")), c(nobs, nvar))
})

test_that("an unsupported sparse.type is rejected", {
  expect_error(
    WriteH5AD(make_srt(), paste0(file_temp(), ".h5ad"), sparse.type = "coo_matrix"),
    "not supported"
  )
})

test_that("a csr_matrix .h5ad can be read", {
  srt <- ReadH5AD(fileh5ad_r)
  counts <- GetAssayData(srt, layer = "counts")

  # Seurat only supports dgCMatrix as counts
  expect_true("dgCMatrix" %in% class(counts))
  expect_equal(dim(counts), c(nvar, nobs))
  expect_equal(counts[true_val_i, true_val_j], true_val)
  expect_equal(rownames(srt), var_names)
  expect_equal(colnames(srt), obs_names)
})

test_that("a csc_matrix .h5ad can be read", {
  srt <- ReadH5AD(fileh5ad_c)
  counts <- GetAssayData(srt, layer = "counts")

  expect_true("dgCMatrix" %in% class(counts))
  expect_equal(dim(counts), c(nvar, nobs))
  expect_equal(counts[true_val_i, true_val_j], true_val)
  expect_equal(rownames(srt), var_names)
  expect_equal(colnames(srt), obs_names)
})
