library(Seurat)
library(MuDataSeurat)
library(Matrix)
library(hdf5r)
library(fs)

# Regression tests for metadata (obs) round-trips.
# Seurat v5 (Assay5) is the baseline.

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

test_that("factor columns round-trip when a category is unused", {
  # Only "a" and "c" are observed, but "b" is a declared level.
  # Truncating the category list by the number of observed codes used to
  # relabel code 2 ("c") as "b".
  celltype <- factor(rep(c("a", "c"), length.out = nobs),
                     levels = c("a", "b", "c"))

  srt <- make_srt()
  srt$celltype <- celltype

  file <- paste0(file_temp(), ".h5ad")
  expect_true(WriteH5AD(srt, file))

  srt2 <- ReadH5AD(file)
  expect_equal(as.character(srt2$celltype), as.character(celltype))
  expect_equal(levels(srt2$celltype), levels(celltype))
})

test_that("factor columns with NAs round-trip", {
  celltype <- factor(c("a", NA, "b", NA, rep("a", nobs - 4)),
                     levels = c("a", "b"))

  srt <- make_srt()
  srt$celltype <- celltype

  file <- paste0(file_temp(), ".h5ad")
  expect_true(WriteH5AD(srt, file))

  srt2 <- ReadH5AD(file)
  expect_equal(as.character(srt2$celltype), as.character(celltype))
  expect_equal(unname(is.na(srt2$celltype)), is.na(celltype))
})

test_that("character columns with NAs round-trip instead of becoming 'NaN'", {
  # AnnData has no nullable-string encoding, so these are stored as
  # categorical. They previously came back as the literal string "NaN".
  batch <- c("b1", NA, "b2", rep("b1", nobs - 3))

  srt <- make_srt()
  srt$batch <- batch

  file <- paste0(file_temp(), ".h5ad")
  expect_true(WriteH5AD(srt, file))

  srt2 <- ReadH5AD(file)
  expect_equal(as.character(srt2$batch), batch)
  expect_true(is.na(as.character(srt2$batch)[2]))
  expect_false(any(as.character(srt2$batch) %in% "NaN"))
})

test_that("character columns without NAs round-trip as strings", {
  batch <- rep(c("b1", "b2"), length.out = nobs)

  srt <- make_srt()
  srt$batch <- batch

  file <- paste0(file_temp(), ".h5ad")
  expect_true(WriteH5AD(srt, file))

  srt2 <- ReadH5AD(file)
  expect_equal(as.character(srt2$batch), batch)
})

test_that("columns containing '/' are renamed consistently in column-order", {
  # HDF5 forbids "/" in object names. Renaming the dataset but leaving the
  # original name in "column-order" produced a file that anndata and mudata
  # could not open at all.
  values <- seq_len(nobs) / 2

  srt <- make_srt()
  srt[["Count (cells/ul)"]] <- values

  file <- paste0(file_temp(), ".h5ad")
  expect_warning(expect_true(WriteH5AD(srt, file)), "does not allow '/'")

  h5 <- H5File$new(file, mode = "r")
  on.exit(h5$close_all())
  written <- names(h5[["obs"]])
  order <- h5attributes(h5[["obs"]])[["column-order"]]

  expect_false("Count (cells/ul)" %in% written)
  expect_true("Count (cells_ul)" %in% written)
  # Every declared column must exist as a dataset, or readers cannot open it.
  expect_true(all(order %in% written))

  srt2 <- ReadH5AD(file)
  expect_equal(unname(srt2[["Count (cells_ul)"]][, 1]), values)
})

test_that("renaming a '/' column does not overwrite a colliding column", {
  slashed <- seq_len(nobs) / 2
  existing <- rep(1, nobs)

  srt <- make_srt()
  srt[["Count (cells/ul)"]] <- slashed
  srt[["Count (cells_ul)"]] <- existing

  file <- paste0(file_temp(), ".h5ad")
  expect_warning(expect_true(WriteH5AD(srt, file)), "does not allow '/'")

  srt2 <- ReadH5AD(file)
  # Both columns survive; the renamed one is suffixed rather than clobbering.
  expect_equal(unname(srt2[["Count (cells_ul)"]][, 1]), existing)
  expect_equal(unname(srt2[["Count (cells_ul).1"]][, 1]), slashed)
})

test_that("stored obs columns replace the ones CreateSeuratObject regenerates", {
  # orig.ident/nCount_*/nFeature_* used to be appended alongside the generated
  # ones as orig.ident.1 etc., leaving the stored values in the .1 columns.
  srt <- make_srt()
  srt$orig.ident <- rep(c("A", "B"), length.out = nobs)

  file <- paste0(file_temp(), ".h5ad")
  expect_true(WriteH5AD(srt, file))

  srt2 <- ReadH5AD(file)
  expect_false(any(grepl("\\.1$", colnames(srt2@meta.data))))
  expect_equal(as.character(srt2$orig.ident), as.character(srt$orig.ident))
  expect_equal(unname(srt2$nCount_RNA), unname(srt$nCount_RNA))
})

test_that("ReadH5MU loads global obs metadata", {
  # Global /obs was read and then discarded, so all meta.data was lost.
  batch <- rep(c("b1", "b2"), length.out = nobs)
  orig <- rep(c("A", "B"), length.out = nobs)

  srt <- make_srt()
  srt$orig.ident <- orig
  srt$batch <- batch

  file <- paste0(file_temp(), ".h5mu")
  expect_true(WriteH5MU(srt, file))

  srt2 <- ReadH5MU(file)
  expect_true("batch" %in% colnames(srt2@meta.data))
  expect_equal(unname(as.character(srt2$batch)), batch)
  expect_equal(unname(as.character(srt2$orig.ident)), orig)
  expect_false(any(grepl("\\.1$", colnames(srt2@meta.data))))
})

test_that("integer columns with NAs round-trip via the nullable encoding", {
  # AnnData has no NA-in-an-integer-array encoding, so such a column is stored
  # as a `values`/`mask` pair. Writing one used to recurse until the stack ran
  # out: the NA-carrying vector was handed straight back to write_matrix().
  srt <- make_srt()
  rank <- c(1L, NA, 3L, NA, 5L, NA, 7L, 8L, NA, 10L)
  srt$rank <- rank

  file <- paste0(file_temp(), ".h5ad")
  expect_true(WriteH5AD(srt, file))

  encoding <- local({
    h5 <- H5File$new(file, mode = "r")
    on.exit(h5$close_all())
    list(
      type = h5attr(h5[["obs"]][["rank"]], "encoding-type"),
      dtype = h5[["obs"]][["rank"]][["values"]]$get_type()$to_text()
    )
  })
  expect_equal(encoding$type, "nullable-integer")
  # The masked-out slots are filled before writing; the fill must not widen the
  # column to a float on disk.
  expect_match(encoding$dtype, "I(32|64)", ignore.case = TRUE)

  back <- ReadH5AD(file)
  expect_equal(unname(back$rank), rank)
  expect_equal(unname(is.na(back$rank)), is.na(rank))
})

test_that("logical columns with NAs round-trip via the nullable encoding", {
  srt <- make_srt()
  flag <- c(TRUE, NA, FALSE, TRUE, NA, FALSE, TRUE, FALSE, NA, TRUE)
  srt$flag <- flag

  file <- paste0(file_temp(), ".h5ad")
  expect_true(WriteH5AD(srt, file))

  encoding <- local({
    h5 <- H5File$new(file, mode = "r")
    on.exit(h5$close_all())
    h5attr(h5[["obs"]][["flag"]], "encoding-type")
  })
  expect_equal(encoding, "nullable-boolean")

  back <- ReadH5AD(file)
  expect_equal(unname(as.logical(back$flag)), flag)
})

test_that("logical columns are stored as a two-value boolean", {
  srt <- make_srt()
  srt$flag <- rep(c(TRUE, FALSE), length.out = nobs)

  file <- paste0(file_temp(), ".h5ad")
  expect_true(WriteH5AD(srt, file))

  labels <- local({
    h5 <- H5File$new(file, mode = "r")
    on.exit(h5$close_all())
    h5[["obs"]][["flag"]]$get_type()$get_labels()
  })
  # hdf5r's default logical type carries a third "NA" level, which h5py reads as
  # uint8 rather than bool -- and anndata refuses a non-boolean nullable mask.
  expect_equal(labels, c("FALSE", "TRUE"))
})
