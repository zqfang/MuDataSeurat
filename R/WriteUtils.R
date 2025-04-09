.mudataversion <- "0.1.0"
.anndataversion <- "0.1.0"
.name <- paste0(getPackageName(), ".r")
.version <- as.character(packageVersion(getPackageName()))

#' @import hdf5r
open_h5 <- function(filename) {
  h5p_create <- H5P_FILE_CREATE$new()
  h5p_create$set_userblock(512)
  H5File$new(filename, mode = "w", file_create_pl = h5p_create)
}

#' @import hdf5r
finalize_mudata <- function(h5) {
  h5$create_attr("encoding-type", "mudata", space = H5S$new("scalar"))
  h5$create_attr("encoding-version", .mudataversion, space = H5S$new("scalar"))
  h5$create_attr("encoder", .name, space = H5S$new("scalar"))
  h5$create_attr("encoder-version", .version, space = H5S$new("scalar"))

  filename <- h5$get_filename()
  h5$close_all()
  h5 <- file(filename, "r+b")
  writeChar(paste0("MuData (format-version=", .mudataversion, ";creator=", .name, ";creator-version=", .version, ")"), h5)
  close(h5)
}

#' @import hdf5r
finalize_anndata_internal <- function(h5) {
  h5$create_attr("encoding-type", "anndata", space = H5S$new("scalar"))
  h5$create_attr("encoding-version", .anndataversion, space = H5S$new("scalar"))
  h5$create_attr("encoder", .name, space = H5S$new("scalar"))
  h5$create_attr("encoder-version", .version, space = H5S$new("scalar"))
}

#' @import hdf5r
finalize_anndata <- function(h5, internal = FALSE) {
  if (internal) {
    finalize_anndata_internal(h5)
  }
  filename <- h5$get_filename()
  h5$close_all()
  h5 <- file(filename, "r+b")
  writeChar(paste0("anndata (format-version=", .mudataversion, ";creator=", .name, ";creator-version=", .version, ")"), h5)
  close(h5)
}

write_dataset <- function(parent, key, obj, scalar = FALSE) {
  dtype <- NULL
  space <- NULL
  if (is.character(obj)) {
    dtype <- H5T_STRING$new(type = "c", size = Inf)
    dtype$set_cset("UTF-8")
  }
  if (scalar) {
    space <- H5S$new("scalar")
  }
  parent$create_dataset(key, obj, dtype = dtype, space = space)
}

write_attribute <- function(obj, name, value, scalar = TRUE) {
  dtype <- NULL
  space <- NULL
  if (is.character(value)) {
    dtype <- H5T_STRING$new(type = "c", size = Inf)
    dtype$set_cset("UTF-8")
  }
  if (length(value) == 1 && scalar) {
    space <- H5S$new("scalar")
  }
  obj$create_attr(name, value, dtype = dtype, space = space)
}

write_matrix <- function(parent, key, mat, storage_sparse_type = "csr_matrix") {
  if (is.matrix(mat) || is.vector(mat) || is.array(mat)) {
    hasna <- anyNA(mat)
    if (hasna && is.double(mat)) {
      # FIXME: extend anndata spec to handle double NAs?
      mat[is.na(mat)] <- NaN
      hasna <- FALSE
    }

    if (hasna && is.character(mat)) {
      # FIXME: anndata spec to handle character with NAs ?
      mat[is.na(mat)] <- NaN
      hasna <- FALSE
    }

    if (!hasna) {
      dset <- write_dataset(parent, key, mat)
      write_attribute(dset, "encoding-type", ifelse(is.character(mat), "string-array", "array"))
      write_attribute(dset, "encoding-version", "0.2.0")
    } else {
      grp <- parent$create_group(key)
      # if a charactor vector with NA, this recursion will lead to stack overflow
      # fixed above
      write_matrix(grp, "values", mat)
      write_matrix(grp, "mask", is.na(mat))
      write_attribute(grp, "encoding-type", ifelse(is.logical(mat), "nullable-boolean", "nullable-integer"))
      write_attribute(grp, "encoding-version", "0.1.0")
    }
  } else if (is.factor(mat)) {
    grp <- parent$create_group(key)
    codes <- as.integer(mat)
    codes[is.na(mat)] <- 0L
    write_matrix(grp, "codes", codes - 1L)
    write_matrix(grp, "categories", levels(mat))
    write_attribute(grp, "ordered", is.ordered(mat))
    write_attribute(grp, "encoding-type", "categorical")
    write_attribute(grp, "encoding-version", "0.2.0")
  } else if (is(mat, "dgCMatrix") || is(mat, "dgRMatrix")) {
    grp <- parent$create_group(key)
    mat0 <- mat
    ## dgCMatrix in R (column-oriented sparse), need to transpose to save as csc via hdf5r
    if (is(mat0, "dgCMatrix") && storage_sparse_type == "csc_matrix") mat0 <- Matrix::t(mat)
    ## dgRMatrix in R (Row-oriented sparse), need to transpose to save as csc via hdf5r
    if (is(mat0, "dgRMatrix") && storage_sparse_type == "csr_matrix") mat0 <- Matrix::t(mat)
    write_dataset(grp, "indptr", mat0@p)
    write_dataset(grp, "data", mat0@x)
    write_attribute(grp, "shape", rev(dim(mat)))
    write_attribute(grp, "encoding-version", "0.1.0")
    if (is(mat0, "dgCMatrix")) {
      write_dataset(grp, "indices", mat0@i)
      write_attribute(grp, "encoding-type", storage_sparse_type) # "csr_matrix")
    } else {
      write_dataset(grp, "indices", mat0@j)
      write_attribute(grp, "encoding-type", storage_sparse_type) # "csc_matrix")
    }
  } else {
    stop("Writing matrices of type ", class(mat), " is not implemented: ", key)
  }
}

write_data_frame <- function(parent, key, attr_df) {
  grp <- parent$create_group(key)
  if (!is.data.frame(attr_df)) { # row names only. Creating a data.frame with duplicated row.names is not possible
    attr_df <- data.frame("_index" = attr_df, check.names = FALSE)
    attr_columns <- character()
  } else {
    attr_columns <- colnames(attr_df)
    attr_df["_index"] <- rownames(attr_df)
  }


  for (col in colnames(attr_df)) {
    # Check if the column is of (Date, POSIXct/POSIXt)
    if (inherits(attr_df[[col]], "Date") ||
      inherits(attr_df[[col]], "POSIXct") ||
      inherits(attr_df[[col]], "POSIXt")) {
      message("Column ", col, " is of datetime type.")
      attr_df[[col]] <- as.character(attr_df[[col]])
    }

    # debug: skip column with all values are NA, it's not allow
    if (all(is.na(attr_df[[col]]))) {
      attr_columns <- attr_columns[attr_columns != col] # remove from column-order
      warning("Skip meta.data column: ", col, ", because of all values are NA.")
      next
    }
    # debug: remove "/" in key, it's not allow
    col2 <- gsub("/", "", col)
    write_matrix(grp, col2, attr_df[[col]])
  }

  # Write attributes
  write_attribute(grp, "_index", "_index")
  write_attribute(grp, "encoding-type", "dataframe")
  write_attribute(grp, "encoding-version", "0.2.0")
  if (length(attr_columns) > 0) {
    write_attribute(grp, "column-order", attr_columns, scalar = FALSE)
  } else {
    # When there are no columns, null buffer can't be written to a file.
    grp$create_attr("column-order", dtype = h5types$H5T_NATIVE_DOUBLE, space = H5S$new("simple", 0, 0))
  }
}

# Only write _index (obs_names or var_names)
write_names <- function(attr_group, attr_names) {
  stype <- H5T_STRING$new(type = "c", size = Inf)
  stype$set_cset("UTF-8")
  attr_group$create_dataset("_index", attr_names, dtype = stype)

  # Write attributes
  attr_group$create_attr("_index", "_index", space = H5S$new("scalar"), dtype = stype)
  attr_group$create_attr("encoding-type", "dataframe", space = H5S$new("scalar"), dtype = stype)
  attr_group$create_attr("encoding-version", "0.2.0", space = H5S$new("scalar"), dtype = stype)
  # When there are no columns, null buffer can't be written to a file.
  attr_group$create_attr("column-order", dtype = h5types$H5T_NATIVE_DOUBLE, space = H5S$new("simple", 0, 0))
}



slot_writer <- function(h5group, mx, name) {
  xt <- mx
  if ("i" %in% slotNames(mx)) {
    sparse_type <- ifelse(class(mx) == "dgCMatrix", "csc_matrix", "csr_matrix")
    # sparse matrix
    if (sparse_type == "csc_matrix") {
      xt <- Matrix::t(mx)
    } # transpose for anndata
    mx_group <- h5group$create_group(name)
    write_sparse_matrix(mx_group, xt, sparse_type)
  } else {
    # dense matrix, create a dataset
    # h5group$create_dataset(name, mx)
    # h5group in the root.
    write_dense_matrix(h5group, xt, name)
  }
}


write_sparse_matrix <- function(root, x, sparse_type) {
  stype <- H5T_STRING$new(type = "c", size = Inf)
  stype$set_cset("UTF-8")
  root$create_dataset("indices", x@i)
  root$create_dataset("indptr", x@p)
  root$create_dataset("data", x@x)
  h5attr(root, "shape") <- dim(x)
  root$create_attr("encoding-type", sparse_type, space = H5S$new("scalar"), dtype = stype)
  root$create_attr("encoding-version", "0.1.0", space = H5S$new("scalar"), dtype = stype)
}

write_dense_matrix <- function(root, x, name) {
  stype <- H5T_STRING$new(type = "c", size = Inf)
  stype$set_cset("UTF-8")
  dense <- root$create_dataset(name, x)
  # h5attr(dense, "shape") <- dim(x)
  dense$create_attr("encoding-type", "array", space = H5S$new("scalar"), dtype = stype)
  dense$create_attr("encoding-version", "0.2.0", space = H5S$new("scalar"), dtype = stype)
}


reshape_scaled_data <- function(mat, var.meta, mat_name = "scale.data") {
  # If only a subset of features was used,
  # this has to be accounted for
  all_mat <- mat
  if (nrow(mat) < nrow(var.meta)) {
    warning(paste0(
      "data values for `", mat_name, "` are computed only for a some features (HVGs).",
      " For it, an array with full var dimension will be recorded as it has to be match the var dimension of the data/counts."
    ))
    all_mat <- matrix(
      ncol = ncol(mat),
      nrow = nrow(var.meta)
    )
    rownames(all_mat) <- rownames(var.meta)
    all_mat[rownames(mat), ] <- mat
  }
  ## don't transpose the dense matrix for anndata (will do the transpose implicity when using hdf5r to write)
  return(all_mat)
}
