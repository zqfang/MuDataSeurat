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
