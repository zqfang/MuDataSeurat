library(Seurat)
library(SeuratObject)
library(MuDataSeurat)
library(Matrix)
library(fs) # for file_temp()

# ConvertToSeuratBPCells()/ConvertToSeuratInMemory() move an object's layers
# between memory and BPCells' on-disk format, leaving everything else about the
# object alone.

nvar <- 20
nobs <- 40

var_names <- paste("var", seq_len(nvar), sep = "-")
obs_names <- paste("obs", seq_len(nobs), sep = "-")

skip_if_no_bpcells <- function() {
  skip_if_not_installed("BPCells")
}

make_counts <- function(seed = 1) {
  set.seed(seed)
  x <- Matrix(matrix(rpois(nvar * nobs, 3), nrow = nvar), sparse = TRUE)
  x <- as(x, "dgCMatrix")
  rownames(x) <- var_names
  colnames(x) <- obs_names
  x
}

# An object with all three layers, feature metadata, variable features and a
# reduction, so the conversions can be checked not to disturb any of it.
make_srt <- function() {
  srt <- CreateSeuratObject(make_counts())
  srt <- NormalizeData(srt, verbose = FALSE)
  srt <- suppressWarnings(FindVariableFeatures(srt, nfeatures = 5, verbose = FALSE))
  srt <- ScaleData(srt, verbose = FALSE)
  assay <- srt[["RNA"]]
  assay[[]] <- data.frame(
    gene_ids = paste0("ENSG", seq_len(nvar)), row.names = var_names
  )
  srt[["RNA"]] <- assay
  srt$batch <- factor(rep(c("p", "q"), length.out = nobs))
  srt[["pca"]] <- CreateDimReducObject(
    embeddings = matrix(
      seq_len(nobs * 2) / 10, nrow = nobs,
      dimnames = list(obs_names, c("pca_1", "pca_2"))
    ),
    key = "pca_", assay = "RNA"
  )
  srt
}

test_that("ConvertToSeuratBPCells moves layers to disk, values preserved", {
  skip_if_no_bpcells()

  srt <- make_srt()
  backed <- ConvertToSeuratBPCells(srt, dir = file_temp())

  expect_true(inherits(LayerData(backed, "counts"), "IterableMatrix"))
  expect_true(inherits(LayerData(backed, "data"), "IterableMatrix"))
  expect_equal(
    unname(as.matrix(LayerData(backed, "counts"))),
    unname(as.matrix(LayerData(srt, "counts")))
  )
  expect_equal(
    unname(as.matrix(LayerData(backed, "data"))),
    unname(as.matrix(LayerData(srt, "data")))
  )
})

test_that("ConvertToSeuratBPCells leaves scale.data in memory by default", {
  skip_if_no_bpcells()

  srt <- make_srt()
  dir <- file_temp()
  backed <- ConvertToSeuratBPCells(srt, dir = dir)

  # scale.data is dense and covers only the variable features, so it gains
  # little from being moved.
  expect_false(inherits(LayerData(backed, "scale.data"), "IterableMatrix"))
  expect_setequal(list.files(file.path(dir, "RNA")), c("counts", "data"))
})

test_that("ConvertToSeuratBPCells converts scale.data when named explicitly", {
  skip_if_no_bpcells()

  srt <- make_srt()
  backed <- ConvertToSeuratBPCells(srt, dir = file_temp(),
                                   layers = "scale.data")

  expect_true(inherits(LayerData(backed, "scale.data"), "IterableMatrix"))
  # The layers not named are left alone.
  expect_false(inherits(LayerData(backed, "counts"), "IterableMatrix"))
})

test_that("ConvertToSeuratBPCells keeps everything that is not a matrix", {
  skip_if_no_bpcells()

  srt <- make_srt()
  backed <- ConvertToSeuratBPCells(srt, dir = file_temp())

  expect_equal(backed[["RNA"]][[]]$gene_ids, srt[["RNA"]][[]]$gene_ids)
  expect_equal(VariableFeatures(backed), VariableFeatures(srt))
  expect_equal(levels(backed$batch), levels(srt$batch))
  expect_equal(Embeddings(backed[["pca"]]), Embeddings(srt[["pca"]]))
  expect_equal(rownames(backed), rownames(srt))
  expect_equal(colnames(backed), colnames(srt))
})

test_that("type casts only the layers it names", {
  skip_if_no_bpcells()

  srt <- make_srt()
  backed <- ConvertToSeuratBPCells(srt, dir = file_temp(),
                             type = c(counts = "uint32_t"))

  expect_equal(BPCells::matrix_type(LayerData(backed, "counts")), "uint32_t")
  # data must stay floating point, and unchanged.
  expect_equal(
    unname(as.matrix(LayerData(backed, "data"))),
    unname(as.matrix(LayerData(srt, "data")))
  )
})

test_that("casting a non-integer layer to an integer type is refused", {
  skip_if_no_bpcells()

  srt <- make_srt()
  # An unnamed type applies to every layer, which would silently truncate the
  # normalized data layer.
  expect_error(
    ConvertToSeuratBPCells(srt, dir = file_temp(), type = "uint32_t"),
    "would corrupt it"
  )
})

test_that("ConvertToSeuratInMemory brings the layers back exactly", {
  skip_if_no_bpcells()

  srt <- make_srt()
  restored <- ConvertToSeuratInMemory(
    ConvertToSeuratBPCells(srt, dir = file_temp())
  )

  for (layer in c("counts", "data")) {
    expect_s4_class(LayerData(restored, layer), "dgCMatrix")
    expect_equal(
      unname(as.matrix(LayerData(restored, layer))),
      unname(as.matrix(LayerData(srt, layer)))
    )
  }
  expect_equal(VariableFeatures(restored), VariableFeatures(srt))
  expect_equal(restored[["RNA"]][[]]$gene_ids, srt[["RNA"]][[]]$gene_ids)
})

test_that("ConvertToSeuratInMemory is a no-op on an in-memory object", {
  srt <- make_srt()
  restored <- ConvertToSeuratInMemory(srt)

  expect_equal(Layers(restored[["RNA"]]), Layers(srt[["RNA"]]))
  expect_equal(
    unname(as.matrix(LayerData(restored, "counts"))),
    unname(as.matrix(LayerData(srt, "counts")))
  )
})

test_that("a converted object survives saveRDS and the session ending", {
  skip_if_no_bpcells()

  srt <- make_srt()
  backed <- ConvertToSeuratBPCells(srt, dir = file_temp())

  rds <- paste0(file_temp(), ".rds")
  saveRDS(backed, rds)
  reloaded <- readRDS(rds)

  expect_true(inherits(LayerData(reloaded, "counts"), "IterableMatrix"))
  expect_equal(
    unname(as.matrix(LayerData(reloaded, "counts"))),
    unname(as.matrix(LayerData(srt, "counts")))
  )
})

test_that("ConvertToSeuratBPCells consolidates a zero-copy object", {
  skip_if_no_bpcells()

  counts <- make_counts()
  file <- paste0(file_temp(), ".h5ad")
  WriteH5AD(CreateSeuratObject(counts), file)

  # Matrices still read from the .h5ad at this point.
  srt <- ReadH5AD(file, backend = "bpcells")
  consolidated <- ConvertToSeuratBPCells(srt, dir = file_temp())

  unlink(file)
  expect_equal(
    unname(as.matrix(LayerData(consolidated, "counts"))),
    unname(as.matrix(counts))
  )
})

test_that("a converted object can still be written back to .h5ad", {
  skip_if_no_bpcells()

  srt <- make_srt()
  backed <- ConvertToSeuratBPCells(srt, dir = file_temp())

  file <- paste0(file_temp(), ".h5ad")
  expect_true(WriteH5AD(backed, file))
  expect_equal(
    unname(as.matrix(LayerData(ReadH5AD(file), "counts"))),
    unname(as.matrix(LayerData(srt, "counts")))
  )
})

test_that("every assay of a multimodal object is converted", {
  skip_if_no_bpcells()

  adt <- make_counts(2)[seq_len(4), , drop = FALSE]
  rownames(adt) <- paste("adt", seq_len(4), sep = "-")
  srt <- CreateSeuratObject(make_counts(), assay = "RNA")
  srt[["ADT"]] <- CreateAssay5Object(counts = adt)

  dir <- file_temp()
  backed <- ConvertToSeuratBPCells(srt, dir = dir)

  expect_setequal(list.files(dir), c("RNA", "ADT"))
  for (assay in c("RNA", "ADT")) {
    expect_true(inherits(LayerData(backed[[assay]], "counts"), "IterableMatrix"))
  }
  expect_equal(
    unname(as.matrix(LayerData(backed[["ADT"]], "counts"))),
    unname(as.matrix(adt))
  )

  # A single assay can be singled out.
  one <- ConvertToSeuratBPCells(srt, dir = file_temp(), assays = "ADT")
  expect_false(inherits(LayerData(one[["RNA"]], "counts"), "IterableMatrix"))
  expect_true(inherits(LayerData(one[["ADT"]], "counts"), "IterableMatrix"))
})

test_that("ConvertToSeuratBPCells refuses a v3 assay", {
  skip_if_no_bpcells()

  srt <- CreateSeuratObject(CreateAssayObject(counts = make_counts()))
  expect_s4_class(srt[["RNA"]], "Assay")
  expect_error(ConvertToSeuratBPCells(srt, dir = file_temp()), "v3 assay")
})

test_that("unknown assays and layers are reported", {
  skip_if_no_bpcells()

  srt <- make_srt()
  expect_error(ConvertToSeuratBPCells(srt, dir = file_temp(), assays = "nope"),
               "Assay not found")
  expect_error(ConvertToSeuratBPCells(srt, dir = file_temp(), layers = "nope"),
               "Layer not found")
  expect_error(ConvertToSeuratInMemory(srt, assays = "nope"), "Assay not found")
})

test_that("an existing matrix directory is not silently overwritten", {
  skip_if_no_bpcells()

  srt <- make_srt()
  dir <- file_temp()
  ConvertToSeuratBPCells(srt, dir = dir)

  expect_error(ConvertToSeuratBPCells(srt, dir = dir))
  expect_silent(
    invisible(ConvertToSeuratBPCells(srt, dir = dir, overwrite = TRUE))
  )
})

test_that("v5 split layers convert to one directory each", {
  skip_if_no_bpcells()

  srt <- CreateSeuratObject(make_counts())
  srt$sample <- rep(c("s1", "s2"), each = nobs / 2)
  srt[["RNA"]] <- split(srt[["RNA"]], f = srt$sample)
  expect_setequal(Layers(srt[["RNA"]]), c("counts.s1", "counts.s2"))

  dir <- file_temp()
  backed <- ConvertToSeuratBPCells(srt, dir = dir)

  expect_setequal(list.files(file.path(dir, "RNA")), c("counts.s1", "counts.s2"))
  for (layer in Layers(srt[["RNA"]])) {
    expect_true(inherits(LayerData(backed, layer), "IterableMatrix"))
    expect_equal(
      unname(as.matrix(LayerData(backed, layer))),
      unname(as.matrix(LayerData(srt, layer)))
    )
  }

  # And the usual v5 workflow still applies once they are back in memory.
  joined <- JoinLayers(ConvertToSeuratInMemory(backed))
  expect_equal(Layers(joined[["RNA"]]), "counts")
})
