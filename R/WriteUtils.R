.mudataversion <- "0.1.0"
.anndataversion <- "0.1.0"
.name <- paste0(getPackageName(), ".r")
.version <- as.character(packageVersion(getPackageName()))

#' @import hdf5r
open_h5 <- function(filename) {
    h5p_create <- H5P_FILE_CREATE$new()
    h5p_create$set_userblock(512)
    H5File$new(filename, mode="w", file_create_pl=h5p_create)
}

#' @import hdf5r
finalize_mudata <- function(h5) {
    h5$create_attr("encoding-type", "mudata", space=H5S$new("scalar"))
    h5$create_attr("encoding-version", .mudataversion, space=H5S$new("scalar"))
    h5$create_attr("encoder", .name, space=H5S$new("scalar"))
    h5$create_attr("encoder-version", .version, space=H5S$new("scalar"))

    filename <- h5$get_filename()
    h5$close_all()
    h5 <- file(filename, "r+b")
    writeChar(paste0("MuData (format-version=", .mudataversion, ";creator=", .name, ";creator-version=", .version, ")"), h5)
    close(h5)
}

#' @import hdf5r
finalize_anndata_internal <- function(h5) {
    h5$create_attr("encoding-type", "anndata", space=H5S$new("scalar"))
    h5$create_attr("encoding-version", .anndataversion, space=H5S$new("scalar"))
    h5$create_attr("encoder", .name, space=H5S$new("scalar"))
    h5$create_attr("encoder-version", .version, space=H5S$new("scalar"))
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


write_data_frame <- function(attr_group, attr_df) {
  attr_columns <- colnames(attr_df)

  attr_df["_index"] <- rownames(attr_df)

  attr_df <- attr_df[,c("_index", attr_columns),drop=FALSE]

  categories <- list()
  # varlenunicode type (hdf5-type)
  sdtype <- H5T_STRING$new(type="c", size=Inf)
  sdtype$set_cset("UTF-8")
  # write dataframe columns
  for (col in colnames(attr_df)) {
    v <- attr_df[[col]]
    if (is.factor(v)) {
      # Write a factor
      categories[[col]] <- levels(v)
      codes <- as.integer(v) - 1
      codes[is.na(codes)] <- -1
      ## categorical array must stored as group
      ## must contain metadata
      #cat_attr = attr_group$create_dataset(col, codes, dtype = h5types$H5T_NATIVE_INT)
      cat_group = attr_group$create_group(col)
      cat_group$create_attr("encoding-type", 'categorical', space = H5S$new("scalar"), dtype=sdtype)
      cat_group$create_attr("encoding-version", '0.2.0', space = H5S$new("scalar"), dtype=sdtype)
      cat_group$create_attr("ordered", FALSE, space = H5S$new("scalar"))
      cat_codes = cat_group$create_dataset("codes", codes, dtype = h5types$H5T_NATIVE_INT)
      cat_codes$create_attr("encoding-type", "array", space=H5S$new("scalar"), dtype=sdtype)
      cat_codes$create_attr("encoding-version", "0.2.0", space=H5S$new("scalar"), dtype=sdtype) 
      cat_items = cat_group$create_dataset("categories", categories[[col]], dtype=sdtype)
      cat_items$create_attr("encoding-type", 'string-array', space = H5S$new("scalar"), dtype=sdtype)
      cat_items$create_attr("encoding-version", '0.2.0', space = H5S$new("scalar"), dtype=sdtype)  
    } else {
      dtype <- NULL
      enc_type = "array"
      if (is.character(v)) {
          enc_type = "string-array"
          dtype = sdtype
      }
      col_attr = attr_group$create_dataset(col, v, dtype=dtype)
      col_attr$create_attr("encoding-type", enc_type, space = H5S$new("scalar"), dtype=sdtype)
      col_attr$create_attr("encoding-version", '0.2.0', space = H5S$new("scalar"), dtype=sdtype)
      # if (col == "_index")
      # {
      #   ## this attribute used by anndata-rs, but it's optional
      #   # index_type: list, range, intervals
      #   col_attr$create_attr("index_type", "list", space = H5S$new("scalar"), dtype=sdtype)
      # }
    }
  }
  # Write attributes
  attr_group$create_attr("_index", "_index", space = H5S$new("scalar"), dtype=sdtype)
  attr_group$create_attr("encoding-type", "dataframe", space = H5S$new("scalar"), dtype=sdtype)
  attr_group$create_attr("encoding-version", "0.2.0", space = H5S$new("scalar"), dtype=sdtype)
  if (length(attr_columns) > 0) {
    attr_group$create_attr("column-order", attr_columns, dtype=sdtype)
  } else {
    # When there are no columns, null buffer can't be written to a file.
    attr_group$create_attr("column-order", dtype=h5types$H5T_NATIVE_DOUBLE, space=H5S$new("simple", 0, 0))
  }
}

# Only write _index (obs_names or var_names)
write_names <- function(attr_group, attr_names) {
  stype <- H5T_STRING$new(type="c", size=Inf)
  stype$set_cset("UTF-8")
  attr_group$create_dataset("_index", attr_names, dtype=stype)

  # Write attributes
  attr_group$create_attr("_index", "_index", space = H5S$new("scalar"), dtype=stype)
  attr_group$create_attr("encoding-type", "dataframe", space = H5S$new("scalar"), dtype=stype)
  attr_group$create_attr("encoding-version", "0.2.0", space = H5S$new("scalar"), dtype=stype)
  # When there are no columns, null buffer can't be written to a file.
  attr_group$create_attr("column-order", dtype=h5types$H5T_NATIVE_DOUBLE, space=H5S$new("simple", 0, 0))

}

write_sparse_matrix <- function(root, x, sparse_type) {
  stype <- H5T_STRING$new(type="c", size=Inf)
  stype$set_cset("UTF-8")
  root$create_dataset("indices", x@i)
  root$create_dataset("indptr", x@p)
  root$create_dataset("data", x@x)
  h5attr(root, "shape") <- dim(x)
  root$create_attr("encoding-type", sparse_type, space=H5S$new("scalar"), dtype=stype)
  root$create_attr("encoding-version", "0.1.0", space=H5S$new("scalar"), dtype=stype)
}

write_dense_matrix <- function(root, x, name) {
  stype <- H5T_STRING$new(type="c", size=Inf)
  stype$set_cset("UTF-8")
  dense = root$create_dataset(name, x)
  #h5attr(dense, "shape") <- dim(x)
  dense$create_attr("encoding-type", "array", space=H5S$new("scalar"), dtype=stype )
  dense$create_attr("encoding-version", "0.2.0", space=H5S$new("scalar"), dtype=stype) 
}


reshape_scaled_data <- function(mat, var.meta, mat_name="scaled.data") {
    # If only a subset of features was used,
    # this has to be accounted for
    all_mat = mat
    if (nrow(mat) < nrow(var.meta)) {
      warning(paste0("data values for `", mat_name,"` are computed only for a some features (HVGs).",
        " For it, an array with full var dimension will be recorded as it has to be match the var dimension of the data/counts."))
      all_mat <- matrix(
        ncol = ncol(mat),
        nrow = nrow(var.meta)
      )
      rownames(all_mat) <- rownames(var.meta)
      all_mat[rownames(mat),] <- mat
    } 
    ## don't transpose the dense matrix for anndata (will do the transpose implicity when using hdf5r to write)
    return(all_mat)
}