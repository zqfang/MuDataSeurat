library(Seurat)
library(MuDataSeurat)
library(Matrix)
library(hdf5r)
library(fs)

# Regression tests for assay selection and layer handling.
# Seurat v5 (Assay5) is the baseline.

nobs <- 10
nvar <- 20

obs_names <- paste("obs", 1:nobs, sep = "-")
var_names <- paste("var", 1:nvar, sep = "-")

make_counts <- function(vars = var_names) {
  x <- rnbinom(n = nobs * length(vars), prob = .95, size = 10)
  x <- Matrix(matrix(x, ncol = nobs), sparse = TRUE)
  colnames(x) <- obs_names
  rownames(x) <- vars
  x
}

make_multimodal <- function() {
  adt_names <- paste("adt", 1:5, sep = "-")
  srt <- CreateSeuratObject(counts = make_counts(), assay = "RNA")
  srt[["ADT"]] <- CreateAssay5Object(counts = make_counts(adt_names))
  srt
}

test_that("WriteH5AD errors on an assay that does not exist", {
  srt <- make_multimodal()
  file <- paste0(file_temp(), ".h5ad")

  # Must not silently fall back to the first assay.
  expect_error(WriteH5AD(srt, file, assay = "rna"), "not found")
})

test_that("WriteH5AD errors when the assay is ambiguous", {
  srt <- make_multimodal()
  file <- paste0(file_temp(), ".h5ad")

  expect_error(WriteH5AD(srt, file), "has to be provided")
})

test_that("WriteH5AD writes the requested assay, not the first one", {
  srt <- make_multimodal()
  file <- paste0(file_temp(), ".h5ad")

  expect_true(WriteH5AD(srt, file, assay = "ADT"))

  srt2 <- ReadH5AD(file)
  expect_equal(nrow(srt2), 5)
  expect_equal(rownames(srt2), paste("adt", 1:5, sep = "-"))
})

test_that("an assay with only a data layer can be written", {
  # The fallback branch used to always pick `counts`, which is NULL here.
  x <- make_counts()
  # Seurat warns about the missing counts layer when deriving nCount/nFeature.
  srt <- suppressWarnings(CreateSeuratObject(CreateAssay5Object(data = x)))

  file <- paste0(file_temp(), ".h5ad")
  expect_true(WriteH5AD(srt, file))

  srt2 <- ReadH5AD(file)
  expect_equal(dim(srt2), c(nvar, nobs))
})

test_that("a multimodal object round-trips through .h5mu", {
  srt <- make_multimodal()
  file <- paste0(file_temp(), ".h5mu")

  expect_true(WriteH5MU(srt, file))

  srt2 <- ReadH5MU(file)
  expect_setequal(Assays(srt2), c("RNA", "ADT"))
  expect_equal(colnames(srt2), obs_names)
  expect_equal(nrow(srt2[["RNA"]]), nvar)
  expect_equal(nrow(srt2[["ADT"]]), 5)
})

test_that(".h5mu contains every group MuData readers require", {
  # These used to be created lazily (or not at all), which produced files that
  # mudata could not open unless multimodal reductions happened to be present.
  srt <- make_multimodal()
  file <- paste0(file_temp(), ".h5mu")
  expect_true(WriteH5MU(srt, file))

  h5 <- H5File$new(file, mode = "r")
  on.exit(h5$close_all())

  required <- c("mod", "obs", "var", "obsm", "varm",
                "obsp", "varp", "obsmap", "varmap", "uns")
  expect_true(all(required %in% names(h5)))
  for (grp in c("obsm", "varm", "obsp", "varp", "obsmap", "varmap")) {
    expect_equal(h5attr(h5[[grp]], "encoding-type"), "dict")
  }
})

test_that(".h5mu obsmap/varmap index the modalities correctly", {
  srt <- make_multimodal()
  file <- paste0(file_temp(), ".h5mu")
  expect_true(WriteH5MU(srt, file))

  h5 <- H5File$new(file, mode = "r")
  on.exit(h5$close_all())

  # Every cell is present in every assay
  expect_equal(as.integer(h5[["obsmap/RNA"]]$read()), 1:nobs)
  expect_equal(as.integer(h5[["obsmap/ADT"]]$read()), 1:nobs)

  # Global var is RNA (20) followed by ADT (5)
  expect_equal(as.integer(h5[["varmap/RNA"]]$read()), c(1:nvar, rep(0, 5)))
  expect_equal(as.integer(h5[["varmap/ADT"]]$read()), c(rep(0, nvar), 1:5))
})

test_that("a single-modality .h5mu can be read", {
  # names(modalities)[2:length()] used to yield c(NA, "RNA") here.
  srt <- CreateSeuratObject(counts = make_counts(), assay = "RNA")
  file <- paste0(file_temp(), ".h5mu")

  expect_true(WriteH5MU(srt, file))

  srt2 <- ReadH5MU(file)
  expect_equal(Assays(srt2), "RNA")
  expect_equal(colnames(srt2), obs_names)
})
