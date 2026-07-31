# Layers that ConvertToSeuratBPCells() leaves alone by default. scale.data is
# dense, usually holds only the variable features, and Seurat's own ScaleData()
# produces it in memory; moving it to a sparse on-disk format would drop its
# zeros and buy very little.
.memory_only_layers <- "scale.data"

resolve_assays <- function(object, assays) {
  available <- names(object@assays)
  if (is.null(assays)) {
    return(available)
  }
  unknown <- setdiff(assays, available)
  if (length(unknown) > 0) {
    stop("Assay not found: ", paste(unknown, collapse = ", "),
         ". Available assays: ", paste(available, collapse = ", "), ".",
         call. = FALSE)
  }
  assays
}

resolve_layers <- function(assay, assay_name, layers, skip = character()) {
  available <- SeuratObject::Layers(assay)
  if (is.null(layers)) {
    return(setdiff(available, skip))
  }
  unknown <- setdiff(layers, available)
  if (length(unknown) > 0) {
    stop("Layer not found in assay ", assay_name, ": ",
         paste(unknown, collapse = ", "), ". Available layers: ",
         paste(available, collapse = ", "), ".", call. = FALSE)
  }
  layers
}

# Seurat v5 splits layers per sample, so a layer can be called "counts.sample1".
# Anything that would be read as a path separator has to go, or the layer would
# be written somewhere other than the directory it was meant to land in.
layer_dirname <- function(layer) {
  gsub("[/\\\\]", "_", layer)
}

# `type` is either one value for every layer or a vector named by layer, so that
# the usual case -- integer counts alongside a normalized `data` layer that must
# stay floating point -- can be expressed as c(counts = "uint32_t").
resolve_type <- function(type, layer) {
  if (is.null(type)) {
    return(NULL)
  }
  if (is.null(names(type))) {
    if (length(type) != 1) {
      stop("type must be a single value or a vector named by layer.",
           call. = FALSE)
    }
    return(type)
  }
  if (!layer %in% names(type)) {
    return(NULL)
  }
  type[[layer]]
}

.integer_matrix_types <- c("uint32_t", "uint64_t")

# Casting normalized values to an integer type truncates them, and BPCells does
# it silently. Catching it here turns a corrupted `data` layer -- the easy
# mistake, since integer counts really are worth casting -- into a refusal that
# names the layer. Only in-memory matrices are checked; an IterableMatrix would
# have to be streamed in full to find out.
check_lossless_cast <- function(mat, type, layer, assay_name) {
  if (is.null(type) || !type %in% .integer_matrix_types || !is(mat, "dgCMatrix")) {
    return(invisible(NULL))
  }
  if (any(mat@x != trunc(mat@x)) || any(mat@x < 0)) {
    stop("Layer ", layer, " of assay ", assay_name, " holds values that are not ",
         "non-negative integers, so casting it to ", type, " would corrupt it.\n",
         "Cast only the layers that hold raw counts, e.g. ",
         "type = c(counts = \"", type, "\").", call. = FALSE)
  }
  invisible(NULL)
}

#' Move a \code{Seurat} object's matrices onto disk
#'
#' Rewrites the layers of a Seurat v5 object into \pkg{BPCells}' on-disk format
#' and returns an object that reads them from there instead of from memory.
#' Everything else -- metadata, feature metadata, variable features, reductions
#' and graphs -- is carried over untouched.
#'
#' Each layer is written to \code{<dir>/<assay>/<layer>}. Layers that are already
#' disk-backed are rewritten too, which is the way to consolidate an object read
#' with \code{ReadH5AD(..., backend = "bpcells")} -- whose matrices still read
#' from the source \code{.h5ad} -- into a self-contained directory.
#'
#' @param object A \code{Seurat} object whose assays are v5 (\code{Assay5}).
#' @param dir Directory to write the matrices under. Created if absent.
#' @param assays Assays to convert. Defaults to all of them.
#' @param layers Layers to convert. Defaults to every layer except
#'   \code{scale.data}, which is dense, usually covers only the variable
#'   features, and gains little from being moved to disk. Naming it explicitly
#'   converts it anyway.
#' @param type Optional \pkg{BPCells} matrix type to cast to. BPCells' bitpacking
#'   compresses integers far better than floats, so raw counts are worth casting
#'   to \code{"uint32_t"} -- but the same cast would truncate a normalized or
#'   scaled layer, so name the layers it applies to:
#'   \code{type = c(counts = "uint32_t")}. A single unnamed value applies to
#'   every layer. Casting an in-memory layer that does not hold non-negative
#'   integers is refused rather than performed silently.
#' @param overwrite Overwrite a matrix directory that already exists.
#'
#' @return The \code{Seurat} object, with the converted layers disk-backed.
#'
#' @seealso \code{\link{ConvertToSeuratInMemory}} for the other direction.
#'
#' @examples
#' \dontrun{
#' seu <- ConvertToSeuratBPCells(seu, dir = "seu_bpcells",
#'                               type = c(counts = "uint32_t"))
#' saveRDS(seu, "seu.rds")  # the matrices stay in seu_bpcells/
#' }
#'
#' @export
ConvertToSeuratBPCells <- function(object, dir, assays = NULL, layers = NULL,
                                   type = NULL, overwrite = FALSE) {
  if (!requireNamespace("BPCells", quietly = TRUE)) {
    stop("ConvertToSeuratBPCells() requires the BPCells package.\n",
         "Install it with remotes::install_github(\"bnprks/BPCells/r\").",
         call. = FALSE)
  }
  if (missing(dir) || !is.character(dir) || length(dir) != 1) {
    stop("dir must be a single path to write the matrices under.", call. = FALSE)
  }

  for (assay_name in resolve_assays(object, assays)) {
    assay <- object[[assay_name]]
    if (!inherits(assay, "Assay5")) {
      stop("Assay ", assay_name, " is a v3 assay; disk-backed matrices require ",
           "a v5 assay. Convert it first, e.g. with ",
           "object[[\"", assay_name, "\"]] <- as(object[[\"", assay_name,
           "\"]], \"Assay5\").", call. = FALSE)
    }

    for (layer in resolve_layers(assay, assay_name, layers,
                                 skip = .memory_only_layers)) {
      mat <- SeuratObject::LayerData(assay, layer)
      # write_matrix_dir() takes an IterableMatrix or a dgCMatrix; a dense
      # scale.data reaches here only when it was asked for by name.
      if (!inherits(mat, "IterableMatrix") && !is(mat, "dgCMatrix")) {
        mat <- methods::as(mat, "dgCMatrix")
      }
      layer_type <- resolve_type(type, layer)
      if (!is.null(layer_type)) {
        check_lossless_cast(mat, layer_type, layer, assay_name)
        mat <- BPCells::convert_matrix_type(mat, layer_type)
      }

      target <- file.path(dir, assay_name, layer_dirname(layer))
      if (!dir.exists(dirname(target))) {
        dir.create(dirname(target), recursive = TRUE)
      }
      backed <- BPCells::write_matrix_dir(mat, dir = target,
                                          overwrite = overwrite)
      SeuratObject::LayerData(object, layer = layer, assay = assay_name) <- backed
    }
  }

  object
}

#' Bring a disk-backed \code{Seurat} object's matrices into memory
#'
#' The inverse of \code{\link{ConvertToSeuratBPCells}}: reads every disk-backed
#' layer into a sparse in-memory matrix. Layers that are already in memory are
#' left alone, so this is safe to call on any object.
#'
#' Note that this materializes the matrices, so the object has to fit in memory
#' afterwards -- that is the point, but it is worth checking against
#' \code{dim(object)} before calling it on something large.
#'
#' @param object A \code{Seurat} object.
#' @param assays Assays to convert. Defaults to all of them.
#' @param layers Layers to convert. Defaults to all of them.
#'
#' @return The \code{Seurat} object, with the converted layers held in memory as
#'   \code{dgCMatrix} objects.
#'
#' @seealso \code{\link{ConvertToSeuratBPCells}} for the other direction.
#'
#' @examples
#' \dontrun{
#' seu <- ReadH5AD("big.h5ad", backend = "bpcells")
#' seu <- ConvertToSeuratInMemory(seu)
#' }
#'
#' @export
ConvertToSeuratInMemory <- function(object, assays = NULL, layers = NULL) {
  for (assay_name in resolve_assays(object, assays)) {
    assay <- object[[assay_name]]
    if (!inherits(assay, "Assay5")) {
      # A v3 assay cannot hold a disk-backed matrix, so there is nothing to do.
      next
    }

    for (layer in resolve_layers(assay, assay_name, layers)) {
      mat <- SeuratObject::LayerData(assay, layer)
      if (!inherits(mat, "IterableMatrix")) {
        next
      }
      SeuratObject::LayerData(object, layer = layer, assay = assay_name) <-
        methods::as(mat, "dgCMatrix")
    }
  }

  object
}
