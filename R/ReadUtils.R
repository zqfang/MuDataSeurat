# Merge metadata read from the file into the meta.data that CreateSeuratObject
# generated. Columns the file provides replace the generated ones instead of
# being appended next to them: CreateSeuratObject recomputes nCount_*/nFeature_*
# from X, which may hold normalised values rather than raw counts, and defaults
# orig.ident to "SeuratProject", so the stored values are the authoritative ones.
add_meta_data <- function(meta_data, new_meta) {
  if (is.null(new_meta) || ncol(new_meta) == 0) {
    return(meta_data)
  }
  kept <- !colnames(meta_data) %in% colnames(new_meta)
  merged <- cbind.data.frame(meta_data[, kept, drop = FALSE], new_meta)
  colnames(merged) <- make.unique(colnames(merged))
  rownames(merged) <- rownames(meta_data)
  merged
}

# CreateSeuratObject() recomputes nCount_<assay>/nFeature_<assay> from the
# assay's counts layer. That costs a full pass over the matrix plus a second
# sparse allocation for the `x > 0` term, which on large objects dominates the
# time spent reading a file. Skip it when it would be wasted or impossible:
#
#   - The stored metadata already carries both columns. Those are the
#     authoritative values and add_meta_data() overwrites the recomputed ones
#     anyway, so computing them achieves nothing. When they are absent it is not
#     waste: skipping would silently drop them from the resulting object.
#   - The assay has no counts layer, which is the case whenever X maps to data
#     or scale.data. Seurat cannot compute the columns from anything else and
#     emits two warnings on its way to not computing them.
#
# `provided` are the column names of every metadata table that will be merged
# into the object's meta.data.
skip_calcn <- function(assay, assay_name, provided) {
  if (!"counts" %in% SeuratObject::Layers(assay)) {
    return(TRUE)
  }
  all(paste0(c("nCount_", "nFeature_"), assay_name) %in% provided)
}

# Evaluate `expr` with Seurat's nCount/nFeature recomputation turned off.
# Older SeuratObject versions do not consult this option, in which case the
# values are computed as before and only the speedup is lost.
without_calcn <- function(expr) {
  old <- options(Seurat.object.assay.calcn = FALSE)
  on.exit(options(old), add = TRUE)
  force(expr)
}

# subset() copies every matrix an assay holds. When the assay already contains
# exactly the requested cells in the requested order -- the usual case for a
# single modality, or for modalities that share all of their observations --
# that copy is a no-op and can be skipped.
subset_cells <- function(assay, cells) {
  if (identical(colnames(assay), cells)) {
    return(assay)
  }
  subset(assay, cells = cells)
}

#' @importFrom hdf5r is_hdf5 H5File
open_and_check_mudata <- function(filename) {
    if (readChar(filename, 6) != "MuData") {
        if (is_hdf5(filename)) {
            warning("The HDF5 file was not created by MuData tooling, we can't guarantee that everything will work correctly", call.=FALSE)
        } else (
            stop("The file is not an HDF5 file", call.=FALSE)
        )
    }
    H5File$new(filename, mode="r")
}

#' @import hdf5r
open_anndata <- function(filename) {
  # PATH/filename.h5mu/mod/rna => read a single modality from the .h5mu file
  path_fragments <- strsplit(filename, "\\.h5mu")[[1]]
  if (length(path_fragments) == 1) {
    h5 <- H5File$new(filename, mode="r")
  } else {
    h5 <- H5File$new(paste0(path_fragments[1], ".h5mu"), mode="r")
    mod_path <- path_fragments[2]
    if (substr(mod_path, 1, 4) != "/mod") {
      mod_path <- paste0("/mod", mod_path)
    }
    h5 <- h5[[mod_path]]
  }
  h5
}

missing_on_read <- function(loc, desc = "") {
  details <- ""
  if (!is.null(desc) && desc != "") {
    details <- paste0("Seurat does not support ", desc, ".")
  }
  warning(paste0("Missing on read: ", loc, ". ", details), call.=FALSE)
}

read_table_encv1 <- function(dataset, set_index = TRUE) {
  columns <- names(dataset)
  columns <- columns[columns != "__categories"]

  col_list <- lapply(columns, function(name) {
    values <- dataset[[name]]$read()
    values_attr <- tryCatch({
      h5attributes(dataset[[name]])
    }, error = function(e) {
      list()
    })
    if (length(values_attr) > 0) {
      if ("categories" %in% names(values_attr)) {
        # Make factors out of categorical data
        ref <- values_attr$categories
        values_labels <- ref$dereference(obj = NULL)[[1]]
        values <- decode_categorical(as.integer(values), values_labels$read())
      }
    }
    values
  })
  table <- data.frame(Reduce(cbind.data.frame, col_list))
  colnames(table) <- columns
  table
}

# Reconstruct a factor from AnnData categorical storage.
# Codes are 0-based indices into `categories`; a negative code denotes NA.
# The full category list is kept as the factor levels, so categories that are
# not used by any observation are preserved rather than shifting the labels of
# the categories that are used.
decode_categorical <- function(codes, categories) {
  codes <- as.integer(codes)
  codes[!is.na(codes) & (codes < 0L | codes >= length(categories))] <- NA
  factor(categories[codes + 1L], levels = categories)
}

read_column <- function(column, etype, eversion) {
  values <- NULL
  if (identical(etype, "categorical")) {
    if (identical(eversion, "0.2.0")) {
      codes <- column[["codes"]]$read()
      categories <- column[["categories"]]$read()
      values <- decode_categorical(codes, categories)
    } else {
      warning(paste0("Cannot recognise encoding-version ", eversion))
    }
  } else if (identical(etype, "nullable-integer") ||
             identical(etype, "nullable-boolean")) {
    # AnnData's nullable encodings store the data and a missing-value mask as
    # two datasets in a group. Whatever sits under a set mask bit is padding the
    # writer chose, so it has to be restored to NA rather than read as a value.
    values <- column[["values"]]$read()
    values[as.logical(column[["mask"]]$read())] <- NA
  } else {
    values <- column$read()
  }
  values
}

read_table_encv2 <- function(dataset, set_index = TRUE) {
  columns <- names(dataset)

  col_list <- lapply(columns, function(name) {

    col_attr <- tryCatch({
      h5attributes(dataset[[name]])
    }, error = function(e) {
      list("encoding-type" = NULL)
    })

    values <- read_column(dataset[[name]], col_attr$`encoding-type`, col_attr$`encoding-version`)

    values
  })
  table <- data.frame(Reduce(cbind.data.frame, col_list))
  colnames(table) <- columns
  table
}

read_table <- function(dataset, set_index = TRUE) {
  if ("H5Group" %in% class(dataset)) {
    # Table is saved as a group rather than a dataset
    dataset_attr <- tryCatch({
      h5attributes(dataset)
    }, error = function(e) {
      list("_index" = "_index")
    })
    indexcol <- "_index"
    if ("_index" %in% names(dataset_attr)) {
      indexcol <- dataset_attr$`_index`
    }

    encv <- "0.1.0"  # some encoding version by default
    if ("encoding-version" %in% names(dataset_attr)) {
      encv <- dataset_attr$`encoding-version`
    }

    if (encv == "0.1.0") {
      table <- read_table_encv1(dataset, set_index)
    } else if (encv == "0.2.0") {
      table <- read_table_encv2(dataset, set_index)
    } else {
      stop(paste0("Encoding version ", encv, " is not recognised."))
    }

    columns <- colnames(table)

    if ((indexcol %in% colnames(table)) && set_index) {
      rownames(table) <- table[,indexcol,drop=TRUE]
      table <- table[,!colnames(table) %in% c(indexcol),drop=FALSE]
    }

    # Fix column order
    if ("column-order" %in% names(dataset_attr)) {
      ordered_columns <- dataset_attr[["column-order"]]
      # Do not consider index as a column
      ordered_columns <- ordered_columns[ordered_columns != indexcol]
      table <- table[,ordered_columns[ordered_columns %in% columns],drop=FALSE]
    }
  } else {
    table <- dataset$read()
    dataset_attr <- h5attributes(dataset)

    indexcol <- "_index"
    if ("_index" %in% names(dataset_attr)) {
      indexcol <- dataset_attr$`_index`
    }

    if ((indexcol %in% colnames(table)) && set_index) {
      rownames(table) <- table[,indexcol,drop=TRUE]
      table <- table[,!colnames(table) %in% c(indexcol),drop=FALSE]
    }
  }
  table
}

# Build a dgCMatrix straight from the buffers read out of the file, bypassing
# Matrix::sparseMatrix(), which makes an extra pass to sort and validate them.
# Returns NULL when the buffers are not in canonical form (indices sorted and
# unique within each slice); sparseMatrix() repairs those, so callers fall back
# to it rather than failing.
#' @import Matrix
new_dgCMatrix <- function(i, p, x, dims) {
  tryCatch(
    new("dgCMatrix",
        i = as.integer(i), p = as.integer(p), x = as.double(x),
        Dim = as.integer(dims), Dimnames = list(NULL, NULL)),
    error = function(e) NULL
  )
}

#' @import Matrix
read_matrix <- function(dataset) {
  if ("data" %in% names(dataset) && "indices" %in% names(dataset) && "indptr" %in% names(dataset)) {
      i <- dataset[["indices"]]$read()
      p <- dataset[["indptr"]]$read()
      x <- dataset[["data"]]$read()

      rowwise <- FALSE
      if ("encoding-type" %in% h5attr_names(dataset)) {
        rowwise <- h5attr(dataset, "encoding-type") == "csr_matrix"
      }

      if ("shape" %in% h5attr_names(dataset)) {
        X_dims <- h5attr(dataset, "shape")
      } else {
        X_dims <- c(length(p) - 1, max(i) + 1)
        if (rowwise) {
          X_dims <- rev(X_dims)
        }
      }

      # The result is always transposed relative to how AnnData stores the
      # matrix, because AnnData is observations x variables and Seurat is
      # variables x observations.
      if (rowwise) {
        # A CSR matrix of shape (n_obs, n_var) and a CSC matrix of shape
        # (n_var, n_obs) have identical indptr/indices/data buffers, so the
        # transpose is a reinterpretation of what was just read rather than a
        # conversion. Going through sparseMatrix() instead would convert CSR to
        # CSC and then Matrix::t() would convert it straight back.
        X <- new_dgCMatrix(i, p, x, rev(X_dims))
        if (!is.null(X)) {
          return(X)
        }
        X <- Matrix::sparseMatrix(j=i, p=p, x=x, dims=X_dims, index1=FALSE)
      } else {
        # CSC input does need a real transpose, but it can at least be built
        # without sparseMatrix()'s extra pass.
        X <- new_dgCMatrix(i, p, x, X_dims)
        if (is.null(X)) {
          X <- Matrix::sparseMatrix(i=i, p=p, x=x, dims=X_dims, index1=FALSE)
        }
      }

      Matrix::t(X)

    } else {
      dataset$read()
    }
}

# `obs` and `var` are the tables at root/obs and root/var. Callers that have
# already read them pass them in: they are also needed to label obsm/varm/obsp,
# and re-reading /obs is expensive once there are many observations.
#' @import Matrix
read_layers_to_assay <- function(root, modalityname="", obs = NULL, var = NULL) {
  X <- read_matrix(root[['X']])

  if (is.null(var)) {
    var <- read_table(root[['var']])
  }
  if (any(grepl("_", rownames(var)))) {
    example_which <- grep("_", rownames(var))[1]
    example_before <- rownames(var)[example_which]
    rownames(var) <- gsub("_", "-", rownames(var))
    example_after <- rownames(var)[example_which]
    warning(paste0("The var_names from modality ", modalityname, " have been renamed as feature names cannot contain '_'.",
      " E.g. ", example_before, " -> ", example_after, "."))
  }

  if (is.null(obs)) {
    obs <- read_table(root[['obs']])
  }
  # NOTE: obs names must NOT be prefixed with the modality name here.
  # ReadH5MU takes the intersection of obs names across modalities to build a
  # single Seurat object, so per-modality prefixes would make it empty.
  colnames(X) <- rownames(obs)
  rownames(X) <- rownames(var)

  raw <- NULL
  if ("raw" %in% names(root)) {
    raw <- root[['raw']]
    raw.X <- read_matrix(raw[['X']])
    raw.var <- read_table(raw[['var']])
    rownames(raw.X) <- rownames(raw.var)
    colnames(raw.X) <- colnames(X)
    if (nrow(raw.X) != nrow(X)) {
      warning(paste0("Only a subset of mod/", modalityname, "/raw/X is loaded, variables (features) that are not present in mod/", modalityname, "/X are discarded."))
      raw.X <- raw.X[rownames(X),]
    }
  }

  layers <- NULL
  custom_layers <- NULL
  if ("layers" %in% names(root)) {
    layers <- lapply(root[['layers']]$names, function(layer_name) {
      layer <- read_matrix(root[['layers']][[layer_name]])
      rownames(layer) <- rownames(X)
      colnames(layer) <- colnames(X)
      layer
    })
    names(layers) <- root[['layers']]$names
    custom_layers <- names(layers)[!names(layers) %in% c("counts")]
    if (length(custom_layers) > 0) {
      missing_on_read(paste0("some of mod/", modalityname, "/layers"), "custom layers, unless labeled 'counts'")
    }
  }

  # Assumptions:
  #   1. X -> counts
  #   2. raw & X -> data & scale.data
  #   3. layers['counts'] & X -> counts & data
  #   4. layers['counts'], raw, X -> counts, data, scale.data
  #
  # These become v5 layers of the same name. A v5 assay names its layers freely,
  # so the mapping no longer has to route matrices into three fixed slots the way
  # the v3 Assay class required.
  counts_as_layer <- !is.null(layers) && "counts" %in% names(layers)
  if (!counts_as_layer && is.null(raw)) {
    # 1
    assay <- SeuratObject::CreateAssay5Object(counts = X)
  } else {
    if (!is.null(raw)) {
      if (counts_as_layer) {
        # 4
        assay <- SeuratObject::CreateAssay5Object(counts = layers[['counts']])
        SeuratObject::LayerData(assay, "data") <- raw.X
        SeuratObject::LayerData(assay, "scale.data") <- X
      } else {
        # 2
        assay <- SeuratObject::CreateAssay5Object(data = raw.X)
        SeuratObject::LayerData(assay, "scale.data") <- X
      }
    } else {
      # 3
      assay <- SeuratObject::CreateAssay5Object(counts = layers[['counts']])
      SeuratObject::LayerData(assay, "data") <- X
    }
  }

  # A v5 assay keeps feature metadata in @meta.data rather than the v3
  # @meta.features slot, and [[<- matches it to the assay's features by name.
  # Assigning a frame with no columns is rejected, hence the guard.
  if (ncol(var) > 0) {
    assay[[]] <- var
  }

  assay
}

read_attr_m <- function(root, attr_name, dim_names = NULL) {
  if (is.null(dim_names)) {
    attr_df <- read_table(root[[attr_name]])
    dim_names <- rownames(attr_df)
  }
  attrm_name <- paste0(attr_name, "m")

  attrm <- list()
  if (attrm_name %in% names(root)) {
    attrm <- lapply(names(root[[attrm_name]]), function(space) {
      dset <- root[[attrm_name]][[space]]
      if (dset$attr_exists("encoding-type") && h5attr(dset, "encoding-type") == "dataframe") {
        missing_on_read(paste0(root$get_obj_name(), attrm_name, "/", space), "additional metadata dataframes")
        mx <- NULL
      } else {
        mx <- t(read_matrix(dset))
        if (dim(mx)[1] == 1) {
          mx <- t(mx)
        }
        rownames(mx) <- dim_names
      }
      mx
    })

    names(attrm) <- names(root[[attrm_name]])
    attrm <- attrm[!sapply(attrm, is.null)]
  }

  attrm
}

#' @import Seurat methods
read_attr_p <- function(root, attr_name, dim_names = NULL) {
  if (is.null(dim_names)) {
    attr_df <- read_table(root[[attr_name]])
    dim_names <- rownames(attr_df)
  }
  attrp_name <- paste0(attr_name, "p")

  attrp <- list()
  if (attrp_name %in% names(root)) {
    attrp <- lapply(names(root[[attrp_name]]), function(graph) {
      mx <- read_matrix(root[[attrp_name]][[graph]])
      rownames(mx) <- dim_names
      colnames(mx) <- dim_names
      mx
    })

    names(attrp) <- names(root[[attrp_name]])
  }

  attrp
}

# For some common reductions,
# there are conventional names for the loadings slots
OBSM2VARM <- list("X_pca" = "PCs", "X_mofa" = "LFs")
