library(Seurat)
library(MuDataSeurat)
library(Matrix)
library(hdf5r)
library(fs)

# Regression tests for the fast paths added to the read and write code.
# Each of them replaces a slower implementation, so what matters here is that
# the results are unchanged, not that they are quick.

nobs <- 10
nvar <- 20

obs_names <- paste("obs", 1:nobs, sep = "-")
var_names <- paste("var", 1:nvar, sep = "-")

make_srt <- function() {
  x <- rnbinom(n = nobs * nvar, prob = .95, size = 10)
  x <- Matrix(matrix(x, ncol = nobs), sparse = TRUE)
  colnames(x) <- obs_names
  rownames(x) <- var_names
  CreateSeuratObject(counts = x)
}


# --- compression ------------------------------------------------------------

test_that("every compression setting produces an identical object", {
  srt <- make_srt()
  srt$score <- seq_len(nobs) / nobs
  srt$group <- factor(rep(c("a", "b"), length.out = nobs))

  reference <- NULL
  for (compression in list("gzip", "none", 0L, 1L, 9L)) {
    file <- paste0(file_temp(), ".h5ad")
    expect_true(WriteH5AD(srt, file, compression = compression))
    back <- ReadH5AD(file)

    expect_equal(dim(back), c(nvar, nobs))
    expect_equal(GetAssayData(back, layer = "counts"),
                 GetAssayData(srt, layer = "counts"))
    expect_equal(back$score, srt$score)
    expect_equal(back$group, srt$group)

    if (is.null(reference)) {
      reference <- GetAssayData(back, layer = "counts")
    } else {
      expect_equal(GetAssayData(back, layer = "counts"), reference)
    }
  }
})

test_that("compression = \"none\" turns the gzip filter off", {
  srt <- make_srt()
  gz <- paste0(file_temp(), ".h5ad")
  none <- paste0(file_temp(), ".h5ad")
  WriteH5AD(srt, gz, compression = "gzip")
  WriteH5AD(srt, none, compression = "none")

  # A filter can only be applied to a chunked dataset, so an uncompressed
  # write is also a contiguous one.
  h5gz <- H5File$new(gz, mode = "r")
  on.exit(h5gz$close_all(), add = TRUE)
  h5none <- H5File$new(none, mode = "r")
  on.exit(h5none$close_all(), add = TRUE)

  layout <- function(h5) as.integer(h5[["X"]][["data"]]$get_create_plist()$get_layout())
  expect_equal(layout(h5gz), as.integer(h5const$H5D_CHUNKED))
  expect_equal(layout(h5none), as.integer(h5const$H5D_CONTIGUOUS))
})

test_that("compression rejects values it cannot honour", {
  srt <- make_srt()
  file <- paste0(file_temp(), ".h5ad")
  expect_error(WriteH5AD(srt, file, compression = "lzf"), "compression must be")
  expect_error(WriteH5AD(srt, file, compression = 42), "compression must be")
  expect_error(WriteH5AD(srt, file, compression = -1), "compression must be")
})

test_that("compression is honoured by WriteH5MU too", {
  srt <- make_srt()
  file <- paste0(file_temp(), ".h5mu")
  expect_true(WriteH5MU(srt, file, compression = "none"))
  back <- ReadH5MU(file)
  expect_equal(GetAssayData(back, layer = "counts"),
               GetAssayData(srt, layer = "counts"))
})


# --- sparse matrix construction ---------------------------------------------

test_that("csr_matrix and csc_matrix storage read back to the same matrix", {
  srt <- make_srt()
  csr <- paste0(file_temp(), ".h5ad")
  csc <- paste0(file_temp(), ".h5ad")
  WriteH5AD(srt, csr, sparse.type = "csr_matrix")
  WriteH5AD(srt, csc, sparse.type = "csc_matrix")

  from_csr <- GetAssayData(ReadH5AD(csr), layer = "counts")
  from_csc <- GetAssayData(ReadH5AD(csc), layer = "counts")

  expect_s4_class(from_csr, "dgCMatrix")
  expect_s4_class(from_csc, "dgCMatrix")
  expect_equal(from_csr, GetAssayData(srt, layer = "counts"))
  expect_equal(from_csc, from_csr)
})

test_that("new_dgCMatrix declines buffers that are not in canonical form", {
  # Indices out of order within a column: sparseMatrix() repairs this, the
  # direct constructor must refuse it so the caller can fall back.
  expect_null(new_dgCMatrix(i = c(2L, 0L), p = c(0L, 2L), x = c(1, 2), dims = c(3L, 1L)))
  # Duplicated index within a column.
  expect_null(new_dgCMatrix(i = c(1L, 1L), p = c(0L, 2L), x = c(1, 2), dims = c(3L, 1L)))

  ok <- new_dgCMatrix(i = c(0L, 2L), p = c(0L, 2L), x = c(1, 2), dims = c(3L, 1L))
  expect_s4_class(ok, "dgCMatrix")
  expect_equal(as.vector(ok), c(1, 0, 2))
})

test_that("a file with unsorted indices still reads correctly", {
  srt <- make_srt()
  file <- paste0(file_temp(), ".h5ad")
  WriteH5AD(srt, file)
  expected <- GetAssayData(srt, layer = "counts")

  # Reverse the indices within the first row, keeping each value with its own
  # index, so the buffers stop being canonical while the matrix they describe
  # stays the same. The direct constructor must refuse them and the fallback
  # must reconstruct the original matrix.
  h5 <- H5File$new(file, mode = "r+")
  indptr <- h5[["X"]][["indptr"]]$read()
  first <- seq_len(indptr[2])
  if (length(first) > 1) {
    idx <- h5[["X"]][["indices"]]
    dat <- h5[["X"]][["data"]]
    i0 <- idx$read()
    x0 <- dat$read()
    idx[first] <- rev(i0[first])
    dat[first] <- rev(x0[first])
  }
  h5$close_all()

  expect_equal(GetAssayData(ReadH5AD(file), layer = "counts"), expected)
})


# --- nCount / nFeature ------------------------------------------------------

test_that("nCount/nFeature come back when /obs carries them", {
  srt <- make_srt()
  expect_true(all(c("nCount_RNA", "nFeature_RNA") %in% colnames(srt@meta.data)))

  file <- paste0(file_temp(), ".h5ad")
  WriteH5AD(srt, file)
  back <- ReadH5AD(file)

  expect_equal(back$nCount_RNA, srt$nCount_RNA)
  expect_equal(back$nFeature_RNA, srt$nFeature_RNA)
})

test_that("nCount/nFeature are still computed when /obs lacks them", {
  # Skipping the recomputation is only safe when the file supplies the values.
  srt <- make_srt()
  stripped <- srt
  stripped@meta.data <- srt@meta.data[, !colnames(srt@meta.data) %in%
                                        c("nCount_RNA", "nFeature_RNA"), drop = FALSE]

  file <- paste0(file_temp(), ".h5ad")
  WriteH5AD(stripped, file)
  back <- ReadH5AD(file)

  expect_true(all(c("nCount_RNA", "nFeature_RNA") %in% colnames(back@meta.data)))
  expect_equal(back$nCount_RNA, srt$nCount_RNA)
  expect_equal(back$nFeature_RNA, srt$nFeature_RNA)
})

test_that("the calcn option is restored even when object creation fails", {
  before <- getOption("Seurat.object.assay.calcn")
  expect_error(without_calcn(stop("boom")), "boom")
  expect_identical(getOption("Seurat.object.assay.calcn"), before)
})


# --- observation alignment --------------------------------------------------

test_that("subset_cells returns the assay untouched when nothing changes", {
  assay <- make_srt()[["RNA"]]
  expect_identical(subset_cells(assay, colnames(assay)), assay)
})

test_that("subset_cells still subsets when the cells differ", {
  assay <- make_srt()[["RNA"]]
  wanted <- colnames(assay)[1:5]
  # Seurat warns when the layer it copies has more cells than the subset.
  got <- suppressWarnings(subset_cells(assay, wanted))
  expect_equal(colnames(got), wanted)
  expect_equal(ncol(got), 5)
})

