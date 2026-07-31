#' @rdname WriteH5MU
setGeneric("WriteH5MU", function(object, file, scale.data=FALSE, sparse.type="csr_matrix", overwrite = TRUE, compression = "gzip") standardGeneric("WriteH5MU"))

#' @rdname WriteH5AD
setGeneric("WriteH5AD", function(object, file, assay = NULL, scale.data=FALSE, sparse.type="csr_matrix", overwrite = TRUE, compression = "gzip") standardGeneric("WriteH5AD"))

#' A helper function to write a modality (an assay) to an .h5mu file
#'
#' @keywords internal
#'
#' @import hdf5r methods
#' @importFrom Matrix t
WriteH5ADHelper <- function(object, assay, root, sparse.type="csr_matrix", global = FALSE, ds_args = list()) {

  mod_object <- Seurat::GetAssay(object, assay)

  # .obs
  obs_names <- colnames(object)
  # There is no local metadata in Seurat objects
  if (global) {
    obs <- object@meta.data
  } else {
    obs <- data.frame(row.names = obs_names)
  }
  write_data_frame(root, "obs", obs, ds_args = ds_args)

  # .var
  if(inherits(mod_object, "Assay5")) {
    var.features <- mod_object@meta.data$var.features
    var.features <- var.features[!is.na(var.features)]
    var_names <- rownames(mod_object)
    # A v5 assay keeps feature metadata in @meta.data, next to the bookkeeping
    # columns VariableFeatures<- writes there. Those two are re-derived from
    # `highly_variable` on read, so they are dropped rather than written out.
    # Everything else is real feature metadata and was previously discarded.
    meta.features <- mod_object[[]]
    meta.features <- meta.features[
      , !colnames(meta.features) %in% c("var.features", "var.features.rank"),
      drop = FALSE
    ]
  }else{
    # assay v4
    var.features = mod_object@var.features
    meta.features <- mod_object@meta.features
    var_names <- rownames(mod_object@meta.features)

  }

  # Define highly variable features, if any
  if (length(var.features) > 0) {
    meta.features$highly_variable <- rownames(meta.features) %in% var.features
    message(paste0(assay, " Added .var['highly_variable'] with highly variable features to meta.features data"))
  }
  
  
  write_data_frame(root, "var", meta.features, ds_args = ds_args)

  # .X, .layers['counts']
  # Assumptions:
  #   1. counts only, or data only -> X
  #   2. counts & data             -> layers['counts'], X = data
  #
  # scale.data is deliberately never written. ScaleData() runs on the variable
  # features only, while AnnData requires X to span every feature, so exporting
  # it meant padding it out to the full var axis with NaN. That padded array is
  # n_features/n_variable_features times the size of the real data -- commonly
  # 10-15x, and dense -- which made it by far the largest thing a write ever
  # allocated, and it is mostly NaN once written, which is not something a
  # reader can compute on anyway. Scaling is one call to recompute downstream
  # (ScaleData() in Seurat, sc.pp.scale in scanpy), so the data is not lost.
  x_names <- list("counts", "data")

  x <- lapply(x_names, function(x_name) {
    x <- NULL
    # assay v4
    if (inherits(mod_object, "Assay") && x_name %in% slotNames(mod_object)) {
      x <- Seurat::GetAssayData(mod_object, x_name)
      if (nrow(x) == 0 || ncol(x) == 0)
        x <- NULL
    }
    # assay v5
    if (inherits(mod_object, "Assay5") && x_name %in% names(mod_object@layers)) {
      x <- Seurat::GetAssayData(mod_object, layer=x_name, assay=assay)
      if (nrow(x) == 0 || ncol(x) == 0)
        x <- NULL
    }
    x
  })
  names(x) <- unlist(x_names)

  if (!is.null(x[["counts"]]) && !is.null(x[["data"]])) {
    # 2
    layers_group <- root$create_group("layers")
    write_attribute(layers_group, "encoding-type", "dict")
    write_attribute(layers_group, "encoding-version", "0.1.0")
    write_matrix(layers_group, "counts", x[["counts"]], sparse.type, ds_args)
    write_matrix(root, "X", x[["data"]], sparse.type, ds_args)
  } else {
    # 1: exactly one of counts/data is present, write that one as X.
    which_x <- which(!vapply(x, is.null, logical(1)))
    if (length(which_x) == 0) {
      stop(paste0("Assay ", assay, " has no data in counts or data to write."))
    }
    write_matrix(root, "X", x[[which_x[1]]], sparse.type, ds_args)
  }

  uns_group <- root$create_group("uns")
  write_attribute(uns_group, "encoding-type", "dict")
  write_attribute(uns_group, "encoding-version", "0.1.0")

  # reductions -> .obsm
  if ('reductions' %in% slotNames(object)) {
    obsm_group <- root$create_group("obsm")
    write_attribute(obsm_group, "encoding-type", "dict")
    write_attribute(obsm_group, "encoding-version", "0.1.0") 
    varm_group <- root$create_group("varm")
    write_attribute(varm_group, "encoding-type", "dict")
    write_attribute(varm_group, "encoding-version", "0.1.0") 

    for (red_name in names(object@reductions)) {
      red <- object@reductions[[red_name]]
      emb <- t(red@cell.embeddings)
      emb_assay <- red@assay.used
      loadings <- red@feature.loadings

      modality_specific <- FALSE
      # Modality-specific reductions can be identified with all their feature names
      # coming from the @assay.used.
      if (!modality_specific) {
        if (!is.null(loadings) && ncol(loadings) == ncol(red)) {
          if (all(rownames(loadings) %in% var_names)) {
            modality_specific <- TRUE
          }
        }
      }

      # Modality-specific reductions in Seurat objects
      # can also start with modality name by the convention used in this package.
      # Multimodal reductions also have the @assay.used set because this is enforced
      # by the current Seurat package.
      if (!is.null(emb_assay) && emb_assay != "" && emb_assay == assay) {
        # Only count reduction as modality-specific
        # if its name can be found in the reduction name or reduction key.
        # This is required since Seurat does require having an existing modality
        # in assay.used, which complicates loading multimodal embeddings.
        # The latter are currently loaded with the default assay set as assay.used.
        # Seurat strips non-alphanumeric characters when deriving a key from a
        # reduction name (umap_harmony -> umapharmony_), so the name has to be
        # sanitized the same way before comparing it against the key.
        if (grepl(tolower(emb_assay), tolower(red_name), fixed = TRUE) ||
            grepl(tolower(emb_assay), tolower(red@key), fixed = TRUE) ||
            grepl(tolower(red_name), tolower(red@key), fixed = TRUE)  ||
            grepl(tolower(gsub('[^[:alnum:]]', '', red_name)), tolower(red@key), fixed = TRUE)
            ) {
          modality_specific <- TRUE
        }

        # Strip away modality name if the embedding starts with it
        if (emb_assay == substr(red_name, 1, nchar(emb_assay))) {
          red_name <- substr(red_name, nchar(emb_assay) + 1, nchar(red_name))
        }
      }

      if (!modality_specific) {
        warning(paste0("Reduction ", red_name, " (key ", red@key, ", assay.used ",
          emb_assay, ") was not recognised as specific to assay ", assay,
          " and will not be written to .obsm."))
        next
      }

      write_matrix(obsm_group, paste0("X_", red_name), emb, ds_args = ds_args)

      # loadings -> .varm
      if (!is.null(loadings) && ncol(loadings) == ncol(red)) {
        varm_key <- red_name
        if (paste0("X_", red_name) %in% names(OBSM2VARM)) {
          varm_key = OBSM2VARM[[paste0("X_", red_name)]]
        }

        # If only a subset of features was used,
        # this has to be accounted for
        if (nrow(loadings) < nrow(meta.features)) {
          warning(paste0("Loadings for ", red_name, " are computed only for some features.",
            " For it, an array with full var dimension will be recorded as it has to be match the var dimension of the data."))
          all_loadings <- matrix(
            ncol = ncol(loadings),
            nrow = nrow(meta.features)
          )
          rownames(all_loadings) <- rownames(meta.features)
          all_loadings[rownames(loadings),] <- loadings
        } else {
          all_loadings <- loadings
        }

        write_matrix(varm_group, varm_key, t(all_loadings), ds_args = ds_args)
      }

      # stdev -> .uns[...]['variance']
      if (length(red@stdev) > 0) {
        if (!red_name %in% names(uns_group)) {
          uns_red <- uns_group$create_group(red_name)
          write_attribute(uns_red, "encoding-type", "dict")
          write_attribute(uns_red, "encoding-version", "0.1.0")
          write_matrix(uns_red, "variance", red@stdev^2, ds_args = ds_args)
        }
      }
    }
  }

  # graphs -> .obsp
  if ('graphs' %in% slotNames(object)) {
    obsp_group <- root$create_group("obsp")
    write_attribute(obsp_group, "encoding-type", "dict")
    write_attribute(obsp_group, "encoding-version", "0.1.0")  
    for (graph_name in names(object@graphs)) {
      graph <- object@graphs[[graph_name]]
      # Only write the graphs with the correct assay.used
      if ('assay.used' %in% slotNames(graph)) {
        if (length(graph@assay.used) > 0 && graph@assay.used == assay) {
          # Strip away modality name if the graph name starts with it:
          # RNA_distances -> distances
          if (assay == substr(graph_name, 1, nchar(graph@assay.used))) {
            graph_name <- substr(graph_name, nchar(graph@assay.used) + 1, nchar(graph_name))
            # Account for _, which is added by ReadH5AD / ReadH5MU
            if (substr(graph_name, 1, 1) == "_") {
              graph_name <- substr(graph_name, 2, nchar(graph_name))
            }
          }
          write_matrix(obsp_group, graph_name, graph, sparse.type, ds_args)
        }
      }
    }
  }

  finalize_anndata_internal(root)

  TRUE
}

#' Write one assay to .h5ad
#'
#' This function writes the data of one of the assays (modalities) of a \code{Seurat} object into an .h5ad file.
#' The behavior of this function if NAs are present is undefined.
#'
#' The following slots are saved: count matrices (`@counts` and `@data`), `@metadata`, `@reductions`, `@feature.loadings`, `@graphs`.
#'
#' @param object \code{Seurat} object.
#' @param file Path to the .h5ad file.
#' @param assay Assay to write; can be omitted if there is a single assay in the object.
#' @param scale.data Deprecated and ignored; \code{scale.data} is never
#'   written. AnnData requires \code{X} to span every feature, so exporting a
#'   matrix scaled on the variable features only meant padding it out with
#'   \code{NaN} to many times its size. Recompute it after reading instead
#'   (\code{ScaleData()} in Seurat, \code{sc.pp.scale} in scanpy).
#' @param sparse.type String, save as csr_matrix or csc_matrix. Note that
#'   \code{csc_matrix} requires transposing every matrix on the way out, since
#'   Seurat stores them column-oriented; \code{csr_matrix} writes them as-is.
#'   Disk-backed (BPCells) matrices support \code{csr_matrix} only, as they
#'   cannot be transposed without being rewritten in full.
#' @param overwrite Boolean value to indicate if to overwrite the \code{file} if it exists (\code{TRUE} by default).
#' @param compression Compression applied to the HDF5 datasets: \code{"gzip"}
#'   (the default, and what previous versions always did), \code{"none"}, or an
#'   integer gzip level between 0 and 9. Compression, not disk I/O, dominates
#'   the time spent writing large objects, so \code{"none"} is several times
#'   faster in exchange for a considerably larger file.
#'
#' @rdname WriteH5AD
#'
#' @import hdf5r
#'
#' @exportMethod WriteH5AD
setMethod("WriteH5AD", "Seurat", function(object, file, assay = NULL, scale.data=FALSE, sparse.type="csr_matrix", overwrite = TRUE, compression = "gzip") {
  warn_scale_data_deprecated(scale.data)
  if (isFALSE(overwrite) && file.exists(file)) {
    stop(paste0("File ", file, " already exists. Use `overwrite = TRUE` to overwrite it or choose a different file name."))
  }
  if (!sparse.type %in% c("csr_matrix", "csc_matrix")) {
    stop(paste0("sparse.type: ", sparse.type, " not supported. Use `csr_matrix` or `csc_matrix`. "))
  }
  check_not_backing_file(object, file)
  ds_args <- resolve_compression(compression)

  h5 <- open_h5(file)

  # When multiple modalities are present,
  # an assay has to be specified.
  # Do not default to Seurat::DefaultAssay(object)
  # as it is not explicit, is hard to reason about,
  # and does not mean anything for MuData.
  if (is.null(assay)) {
    if (length(object@assays) > 1) {
      h5$close()
      stop(paste0(
        "An assay to be written has to be provided, one of: ",
        paste(names(object@assays), collapse = ", "),
        ".\nUse WriteH5MU() to write all the modalities."
      ))
    }
    assay <- names(object@assays)[1]
  } else if (!assay %in% names(object@assays)) {
    # Never fall back to another assay: writing a different assay than the one
    # that was asked for would silently produce a file with the wrong data.
    h5$close()
    stop(paste0(
      "Assay ", assay, " not found. Available assays: ",
      paste(names(object@assays), collapse = ", "),
      "."
    ))
  }

  # "Global" attributes such as metadata have to be written
  WriteH5ADHelper(object, assay, h5, sparse.type, global = TRUE, ds_args = ds_args)

  finalize_anndata(h5)

  invisible(TRUE)
})

#' Create an .h5mu file with data from a \code{\link{Seurat}} object
#'
#' Save \code{\link{Seurat}} object to .h5mu file.
#' The behavior of this function if NAs are present is undefined.
#'
#' The following slots are saved: count matrices (`@counts` and `@data`), `@metadata`, `@reductions`, `@feature.loadings`, `@graphs`.
#'
#' @param object \code{Seurat} object.
#' @param file Path to the .h5mu file.
#' @param scale.data Deprecated and ignored; \code{scale.data} is never
#'   written. AnnData requires \code{X} to span every feature, so exporting a
#'   matrix scaled on the variable features only meant padding it out with
#'   \code{NaN} to many times its size. Recompute it after reading instead
#'   (\code{ScaleData()} in Seurat, \code{sc.pp.scale} in scanpy).
#' @param sparse.type String, save as csr_matrix or csc_matrix. Note that
#'   \code{csc_matrix} requires transposing every matrix on the way out, since
#'   Seurat stores them column-oriented; \code{csr_matrix} writes them as-is.
#'   Disk-backed (BPCells) matrices support \code{csr_matrix} only, as they
#'   cannot be transposed without being rewritten in full.
#' @param overwrite Boolean value to indicate if to overwrite the \code{file} if it exists (\code{TRUE} by default).
#' @param compression Compression applied to the HDF5 datasets: \code{"gzip"}
#'   (the default, and what previous versions always did), \code{"none"}, or an
#'   integer gzip level between 0 and 9. Compression, not disk I/O, dominates
#'   the time spent writing large objects, so \code{"none"} is several times
#'   faster in exchange for a considerably larger file.
#'
#' @rdname WriteH5MU
#'
#' @import hdf5r methods
#'
#' @exportMethod WriteH5MU
setMethod("WriteH5MU", "Seurat", function(object, file, scale.data=FALSE, sparse.type="csr_matrix", overwrite=TRUE, compression = "gzip") {
  warn_scale_data_deprecated(scale.data)
  if (!sparse.type %in% c("csr_matrix", "csc_matrix")) {
    stop(paste0("sparse.type: ", sparse.type, " not supported. Use `csr_matrix` or `csc_matrix`. "))
  }
  check_not_backing_file(object, file)
  ds_args <- resolve_compression(compression)
  h5 <- open_h5(file)
  # .obs
  obs <- object@meta.data

  write_data_frame(h5, "obs", obs, ds_args = ds_args)

  modalities <- Seurat::Assays(object)

  h5mod <- h5$create_group("mod")
  h5mod$create_attr("mod-order", modalities)
  var_names <- lapply(modalities, function(mod) {
    mod_group <- h5$create_group(paste0("mod/", mod))

    WriteH5ADHelper(object, mod, mod_group, sparse.type, ds_args = ds_args)

    mod_object <- object[[mod]]
    rownames(mod_object)
  })
  names(var_names) <- modalities
  write_data_frame(h5, "var", do.call(c, var_names), ds_args = ds_args)

  write_mod_maps(h5, modalities, nrow(obs), var_names, ds_args = ds_args)

  uns_group <- h5$create_group("uns")
  write_attribute(uns_group, "encoding-type", "dict")
  write_attribute(uns_group, "encoding-version", "0.1.0")

  # obsm/varm/obsp/varp have to exist even when empty: MuData readers expect all
  # of them, and creating them lazily produced files that could not be opened.
  obsm_group <- h5$create_group("obsm")
  write_attribute(obsm_group, "encoding-type", "dict")
  write_attribute(obsm_group, "encoding-version", "0.1.0")
  varm_group <- h5$create_group("varm")
  write_attribute(varm_group, "encoding-type", "dict")
  write_attribute(varm_group, "encoding-version", "0.1.0")
  obsp_group <- h5$create_group("obsp")
  write_attribute(obsp_group, "encoding-type", "dict")
  write_attribute(obsp_group, "encoding-version", "0.1.0")
  varp_group <- h5$create_group("varp")
  write_attribute(varp_group, "encoding-type", "dict")
  write_attribute(varp_group, "encoding-version", "0.1.0")

  # reductions -> .obsm
  # Reductions starting with modality name
  # that corresponds to the assay.used value
  # will be stored in .obsm slots of individual modalities:
  # RNAUMAP -> /mod/RNA/obsm/UMAP
  if ('reductions' %in% slotNames(object)) {
    for (red_name in names(object@reductions)) {
      red <- object@reductions[[red_name]]
      emb <- t(red@cell.embeddings)
      assay_emb <- red@assay.used # assay name which reduction constructed from. 'RNA', 'ADT', 'SCT' etc.
      loadings <- red@feature.loadings
      # reduction.name => red_name , reduction.key => red@key

      modality_specific <- FALSE
      # Modality-specific reductions can be identified with all their feature names
      # coming from the @assay.used.
      if (!modality_specific) {
        if (!is.null(loadings) && ncol(loadings) == ncol(red)) {
          if (all(rownames(loadings) %in% var_names[[assay_emb]])) {
            modality_specific <- TRUE
          }
        }
      }

      # Modality-specific reductions in Seurat objects
      # can also start with modality name by the convention used in this package.
      # Multimodal reductions also have the @assay.used set because this is enforced
      # by the current Seurat package.
      if (!is.null(assay_emb) && assay_emb != "" && assay_emb %in% modalities) {
        # Only count reduction as modality-specific
        # if its name can be found in the reduction name or reduction key.
        # This is required since Seurat does require having an existing modality
        # in assay.used, which complicates loading multimodal embeddings.
        # The latter are currently loaded with the default assay set as assay.used.
        # Seurat strips non-alphanumeric characters when deriving a key from a
        # reduction name (umap_harmony -> umapharmony_), so the name has to be
        # sanitized the same way before comparing it against the key.
        if (grepl(tolower(assay_emb), tolower(red_name), fixed = TRUE) ||
            grepl(tolower(assay_emb), tolower(red@key), fixed = TRUE) ||
            grepl(tolower(red_name), tolower(red@key), fixed = TRUE)  ||
            grepl(tolower(gsub('[^[:alnum:]]', '', red_name)), tolower(red@key), fixed = TRUE)
            ) {
          modality_specific <- TRUE
        }


        # Strip away modality name if the embedding starts with it
        if (assay_emb == substr(red_name, 1, nchar(assay_emb))) {
          red_name <- substr(red_name, nchar(assay_emb) + 1, nchar(red_name))
        }
      }

      if (modality_specific) {
        next
      }

      write_matrix(obsm_group, paste0("X_", red_name), emb, ds_args = ds_args)

      # loadings -> .varm
      if (!is.null(loadings) && ncol(loadings) == ncol(red)) {
        varm_key <- red_name
        if (paste0("X_", red_name) %in% names(OBSM2VARM)) {
          varm_key <- OBSM2VARM[[paste0("X_", red_name)]]
        }

        # If only a subset of features was used,
        # this has to be accounted for
        var_names_for_loadings <- do.call(c, var_names)

        if (nrow(loadings) < length(var_names_for_loadings)) {
          warning(paste0("Loadings for ", red_name, " are computed only for a some features.",
            " For it, an array with full var dimension will be recorded as it has to be match the var dimension of the data."))
          all_loadings <- matrix(
            ncol = ncol(loadings),
            nrow = length(var_names_for_loadings)
          )
          rownames(all_loadings) <- var_names_for_loadings
          all_loadings[rownames(loadings),] <- loadings
        } else {
          all_loadings <- loadings
        }

        write_matrix(varm_group, varm_key, t(all_loadings), ds_args = ds_args)
      }

      # stdev -> .uns[...]['variance']
      # Modality-specific reductions have already been written by WriteH5ADHelper
      # and skipped above, so only multimodal reductions reach this point.
      if (length(red@stdev) > 0) {
        if (!red_name %in% names(uns_group)) {
          uns <- uns_group$create_group(red_name)
          write_attribute(uns, "encoding-type", "dict")
          write_attribute(uns, "encoding-version", "0.1.0")
        } else {
          uns <- uns_group[[red_name]]
        }
        write_matrix(uns, "variance", red@stdev^2, ds_args = ds_args)
      }
    }
  }

  # graphs -> .obsp
  if ('graphs' %in% slotNames(object)) {
    for (graph_name in names(object@graphs)) {
      graph <- object@graphs[[graph_name]]

      # Only write the graphs with no (correct) assay.used
      graph_no_assay <- FALSE
      if (!'assay.used' %in% slotNames(graph)) {
        graph_no_assay <- TRUE
      } else {
        if (length(graph@assay.used) < 1) {
          graph_no_assay <- TRUE
        }
        else if (!graph@assay.used %in% modalities) {
          graph_no_assay <- TRUE
        }
      }

      if (graph_no_assay) {
        write_matrix(obsp_group, graph_name, graph, sparse.type, ds_args)
      }
    }
  }

  finalize_mudata(h5)

  invisible(TRUE)
})
