library(Seurat)
library(SeuratObject)
library(MuDataSeurat)
library(Matrix)
library(hdf5r)
library(fs)

# scale.data is never exported.
#
# ScaleData() runs on the variable features only, but AnnData requires X to span
# every feature, so writing it meant padding the dense scaled matrix out to the
# full var axis with NaN -- many times the size of the real data, and mostly NaN
# once written. Scaling is recomputed in one call downstream instead.

nobs <- 12
nvar <- 10
nhvg <- 4

obs_names <- paste("obs", seq_len(nobs), sep = "-")
var_names <- paste("var", seq_len(nvar), sep = "-")
hvg_names <- var_names[seq_len(nhvg)]

# An object carrying counts, data, and a scale.data covering only the HVGs.
make_scaled_srt <- function() {
  set.seed(11)
  counts <- Matrix(matrix(rpois(nvar * nobs, 5), nrow = nvar), sparse = TRUE)
  counts <- as(counts, "dgCMatrix")
  dimnames(counts) <- list(var_names, obs_names)

  srt <- CreateSeuratObject(counts = counts)
  srt <- SetAssayData(srt, layer = "data", new.data = counts)
  scaled <- matrix(rnorm(nhvg * nobs), nrow = nhvg,
                   dimnames = list(hvg_names, obs_names))
  srt <- SetAssayData(srt, layer = "scale.data", new.data = scaled)
  VariableFeatures(srt) <- hvg_names
  srt
}

h5_names <- function(file, path = NULL) {
  h5 <- H5File$new(file, mode = "r")
  on.exit(h5$close_all())
  if (is.null(path)) names(h5) else names(h5[[path]])
}


test_that("scale.data is not written to .h5ad", {
  srt <- make_scaled_srt()
  expect_true("scale.data" %in% Layers(srt[["RNA"]]))

  file <- paste0(file_temp(), ".h5ad")
  expect_true(WriteH5AD(srt, file))

  # counts and data are both present, so counts goes to layers and data is X.
  expect_setequal(h5_names(file, "layers"), "counts")
  expect_false("raw" %in% h5_names(file))
})

test_that("X is the data layer, not the scaled one", {
  srt <- make_scaled_srt()
  file <- paste0(file_temp(), ".h5ad")
  WriteH5AD(srt, file)

  back <- ReadH5AD(file)
  expect_false("scale.data" %in% Layers(back[["RNA"]]))
  expect_setequal(Layers(back[["RNA"]]), c("counts", "data"))
  expect_equal(LayerData(back[["RNA"]], "counts"),
               LayerData(srt[["RNA"]], "counts"))
  expect_equal(LayerData(back[["RNA"]], "data"),
               LayerData(srt[["RNA"]], "data"))
})

test_that("X spans every feature and holds no NaN", {
  # The padded export used to make X mostly NaN; nothing should be padded now.
  srt <- make_scaled_srt()
  file <- paste0(file_temp(), ".h5ad")
  WriteH5AD(srt, file)

  h5 <- H5File$new(file, mode = "r")
  on.exit(h5$close_all())
  expect_equal(as.integer(h5attr(h5[["X"]], "shape")), c(nobs, nvar))
  expect_false(any(is.nan(h5[["X"]][["data"]]$read())))
})

test_that("scale.data is not written to .h5mu either", {
  srt <- make_scaled_srt()
  file <- paste0(file_temp(), ".h5mu")
  expect_true(WriteH5MU(srt, file))

  expect_setequal(h5_names(file, "mod/RNA/layers"), "counts")
  back <- ReadH5MU(file)
  expect_false("scale.data" %in% Layers(back[["RNA"]]))
})

test_that("an assay holding only scale.data is refused rather than written", {
  srt <- make_scaled_srt()
  srt[["RNA"]]@layers$counts <- NULL
  srt[["RNA"]]@layers$data <- NULL

  file <- paste0(file_temp(), ".h5ad")
  expect_error(WriteH5AD(srt, file), "no data in counts or data")
})


# --- the deprecated argument -------------------------------------------------

test_that("scale.data = TRUE warns and is still ignored", {
  srt <- make_scaled_srt()
  file <- paste0(file_temp(), ".h5ad")

  expect_warning(WriteH5AD(srt, file, scale.data = TRUE),
                 "deprecated and ignored")
  expect_setequal(h5_names(file, "layers"), "counts")

  mu <- paste0(file_temp(), ".h5mu")
  expect_warning(WriteH5MU(srt, mu, scale.data = TRUE),
                 "deprecated and ignored")
  expect_setequal(h5_names(mu, "mod/RNA/layers"), "counts")
})

test_that("scale.data = FALSE stays silent, since that is what happens anyway", {
  srt <- make_scaled_srt()
  file <- paste0(file_temp(), ".h5ad")
  # Existing scripts and the README pass scale.data = FALSE; it must keep
  # working without complaint.
  expect_no_warning(WriteH5AD(srt, file, scale.data = FALSE))
  expect_no_warning(WriteH5MU(srt, paste0(file_temp(), ".h5mu"),
                              scale.data = FALSE))
})

test_that("scale.data is still accepted positionally", {
  # It keeps its slot in the signature so positional calls do not shift.
  srt <- make_scaled_srt()
  file <- paste0(file_temp(), ".h5ad")
  expect_true(WriteH5AD(srt, file, "RNA", FALSE, "csr_matrix", TRUE))
})


# --- objects that never had scale.data ---------------------------------------

test_that("counts-only objects are unaffected", {
  set.seed(3)
  counts <- Matrix(matrix(rpois(nvar * nobs, 5), nrow = nvar), sparse = TRUE)
  counts <- as(counts, "dgCMatrix")
  dimnames(counts) <- list(var_names, obs_names)
  srt <- CreateSeuratObject(counts = counts)

  file <- paste0(file_temp(), ".h5ad")
  WriteH5AD(srt, file)
  expect_false("layers" %in% h5_names(file))

  back <- ReadH5AD(file)
  expect_equal(LayerData(back[["RNA"]], "counts"), counts)
})
