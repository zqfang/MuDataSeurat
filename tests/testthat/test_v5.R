library(Seurat)
library(SeuratObject)
library(MuDataSeurat)
library(Matrix)
library(hdf5r)
library(fs) # for file_temp()

# Reading always builds a Seurat v5 assay. The X/raw/layers cases enumerated in
# read_layers_to_assay() map onto v5 layers of the same name, and feature
# metadata lives in the assay's @meta.data rather than the v3 @meta.features.

nobs <- 8
nvar <- 5

obs_names <- paste("obs", seq_len(nobs), sep = "-")
var_names <- paste("var", seq_len(nvar), sep = "-")

make_matrix <- function(seed) {
  set.seed(seed)
  x <- Matrix(matrix(rpois(nvar * nobs, 3), nrow = nvar), sparse = TRUE)
  x <- as(x, "dgCMatrix")
  rownames(x) <- var_names
  colnames(x) <- obs_names
  x
}

# Build an .h5ad directly, so the raw/layers combinations that this package's
# own writer never emits can still be exercised on read.
write_h5ad_parts <- function(file, X, layers = NULL, raw = NULL, var = NULL) {
  h5 <- MuDataSeurat:::open_h5(file)
  if (is.null(var)) {
    var <- data.frame(row.names = var_names)
  }
  MuDataSeurat:::write_data_frame(h5, "obs", data.frame(row.names = obs_names))
  MuDataSeurat:::write_data_frame(h5, "var", var)
  MuDataSeurat:::write_matrix(h5, "X", X)
  if (!is.null(layers)) {
    grp <- h5$create_group("layers")
    MuDataSeurat:::write_attribute(grp, "encoding-type", "dict")
    MuDataSeurat:::write_attribute(grp, "encoding-version", "0.1.0")
    for (name in names(layers)) {
      MuDataSeurat:::write_matrix(grp, name, layers[[name]])
    }
  }
  if (!is.null(raw)) {
    grp <- h5$create_group("raw")
    MuDataSeurat:::write_matrix(grp, "X", raw)
    MuDataSeurat:::write_data_frame(
      grp, "var", data.frame(row.names = var_names)
    )
  }
  MuDataSeurat:::finalize_anndata(h5, internal = TRUE)
  file
}

# Layers compare by value; dimnames are asserted separately where they matter.
expect_layer <- function(object, layer, expected) {
  expect_equal(
    unname(as.matrix(LayerData(object, layer))),
    unname(as.matrix(expected))
  )
}

test_that("ReadH5AD builds a v5 assay", {
  file <- paste0(file_temp(), ".h5ad")
  write_h5ad_parts(file, make_matrix(1))

  srt <- ReadH5AD(file)
  expect_s4_class(srt[["RNA"]], "Assay5")
  expect_false(inherits(srt[["RNA"]], "Assay"))
})

test_that("X alone becomes a counts layer", {
  file <- paste0(file_temp(), ".h5ad")
  counts <- make_matrix(1)
  write_h5ad_parts(file, counts)

  srt <- ReadH5AD(file)
  expect_equal(Layers(srt[["RNA"]]), "counts")
  expect_layer(srt, "counts", counts)
})

test_that("layers['counts'] plus X become counts and data layers", {
  file <- paste0(file_temp(), ".h5ad")
  counts <- make_matrix(1)
  data <- make_matrix(2)
  write_h5ad_parts(file, data, layers = list(counts = counts))

  srt <- ReadH5AD(file)
  expect_setequal(Layers(srt[["RNA"]]), c("counts", "data"))
  expect_layer(srt, "counts", counts)
  expect_layer(srt, "data", data)
})

test_that("raw plus X become data and scale.data layers", {
  file <- paste0(file_temp(), ".h5ad")
  scaled <- make_matrix(1)
  raw <- make_matrix(2)
  write_h5ad_parts(file, scaled, raw = raw)

  srt <- ReadH5AD(file)
  expect_setequal(Layers(srt[["RNA"]]), c("data", "scale.data"))
  expect_layer(srt, "data", raw)
  expect_layer(srt, "scale.data", scaled)
})

test_that("layers['counts'], raw and X become all three layers", {
  file <- paste0(file_temp(), ".h5ad")
  counts <- make_matrix(1)
  raw <- make_matrix(2)
  scaled <- make_matrix(3)
  write_h5ad_parts(file, scaled, layers = list(counts = counts), raw = raw)

  srt <- ReadH5AD(file)
  expect_setequal(Layers(srt[["RNA"]]), c("counts", "data", "scale.data"))
  expect_layer(srt, "counts", counts)
  expect_layer(srt, "data", raw)
  expect_layer(srt, "scale.data", scaled)
})

test_that("feature metadata lands in the assay's meta.data, unduplicated", {
  file <- paste0(file_temp(), ".h5ad")
  var <- data.frame(
    gene_ids = paste0("ENSG", seq_len(nvar)),
    score = seq_len(nvar) / 10,
    row.names = var_names
  )
  write_h5ad_parts(file, make_matrix(1), var = var)

  srt <- ReadH5AD(file)
  metafeatures <- srt[["RNA"]][[]]
  expect_equal(rownames(metafeatures), var_names)
  # /var used to be attached twice: once by read_layers_to_assay() and again by
  # a cbind in ReadH5AD, which duplicated every column.
  expect_equal(sum(colnames(metafeatures) == "gene_ids"), 1L)
  expect_equal(metafeatures$gene_ids, var$gene_ids)
  expect_equal(metafeatures$score, var$score)
})

test_that("feature metadata survives a write/read round trip", {
  counts <- make_matrix(1)
  srt <- CreateSeuratObject(counts)
  assay <- srt[["RNA"]]
  assay[[]] <- data.frame(
    gene_ids = paste0("ENSG", seq_len(nvar)),
    score = seq_len(nvar) / 10,
    row.names = var_names
  )
  srt[["RNA"]] <- assay
  VariableFeatures(srt) <- var_names[c(1, 3)]

  file <- paste0(file_temp(), ".h5ad")
  WriteH5AD(srt, file)

  # A v5 assay stores feature metadata in @meta.data; the writer used to read it
  # from the v3 slot only and so wrote an empty /var.
  var_columns <- local({
    h5 <- H5File$new(file, mode = "r")
    on.exit(h5$close_all())
    names(h5[["var"]])
  })
  expect_true(all(c("gene_ids", "score") %in% var_columns))

  back <- ReadH5AD(file)
  expect_equal(back[["RNA"]][[]]$gene_ids, paste0("ENSG", seq_len(nvar)))
  expect_equal(back[["RNA"]][[]]$score, seq_len(nvar) / 10)
  expect_equal(VariableFeatures(back[["RNA"]]), var_names[c(1, 3)])
})

test_that("the var.features bookkeeping columns are not written to /var", {
  srt <- CreateSeuratObject(make_matrix(1))
  VariableFeatures(srt) <- var_names[c(2, 4)]

  file <- paste0(file_temp(), ".h5ad")
  WriteH5AD(srt, file)

  h5 <- H5File$new(file, mode = "r")
  on.exit(h5$close_all())
  # VariableFeatures<- writes var.features/var.features.rank into @meta.data.
  # They are re-derived from highly_variable on read, so they stay out of /var.
  expect_false(
    any(c("var.features", "var.features.rank") %in% names(h5[["var"]]))
  )
  expect_true("highly_variable" %in% names(h5[["var"]]))
})

test_that("every modality of an .h5mu is read as a v5 assay", {
  rna <- make_matrix(1)
  adt <- make_matrix(2)[seq_len(3), , drop = FALSE]
  rownames(adt) <- paste("adt", seq_len(3), sep = "-")

  srt <- CreateSeuratObject(rna, assay = "RNA")
  srt[["ADT"]] <- CreateAssay5Object(counts = adt)

  rna_assay <- srt[["RNA"]]
  rna_assay[[]] <- data.frame(
    gene_ids = paste0("ENSG", seq_len(nvar)), row.names = var_names
  )
  srt[["RNA"]] <- rna_assay

  file <- paste0(file_temp(), ".h5mu")
  WriteH5MU(srt, file)

  back <- ReadH5MU(file)
  expect_setequal(names(back@assays), c("RNA", "ADT"))
  expect_s4_class(back[["RNA"]], "Assay5")
  expect_s4_class(back[["ADT"]], "Assay5")
  expect_equal(back[["RNA"]][[]]$gene_ids, paste0("ENSG", seq_len(nvar)))
  expect_layer(back[["ADT"]], "counts", adt)
})
