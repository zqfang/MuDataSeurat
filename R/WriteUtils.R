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
  writeChar(paste0("anndata (format-version=", .anndataversion, ";creator=", .name, ";creator-version=", .version, ")"), h5)
  close(h5)
}

# hdf5r compresses every dataset it creates with gzip by default. On a large
# object that compression, not the I/O, is what dominates the time spent
# writing, so it has to be controllable. Translate the user-facing `compression`
# value into the arguments create_dataset() expects; the result is passed down
# unchanged as `ds_args` so it is parsed once per file rather than per dataset.
#
# Without a filter there is nothing for chunking to buy on these write-once,
# read-in-full datasets, so "no compression" also means a contiguous layout,
# which is the fastest of the three.
resolve_compression <- function(compression) {
  invalid <- function() {
    stop("compression must be \"gzip\", \"none\", or an integer between 0 and 9.",
         call. = FALSE)
  }
  if (is.character(compression)) {
    if (length(compression) != 1) invalid()
    level <- switch(compression, gzip = NA_integer_, none = 0L, invalid())
  } else {
    if (length(compression) != 1) invalid()
    level <- suppressWarnings(as.integer(compression))
    if (is.na(level) || level < 0L || level > 9L) invalid()
  }
  if (is.na(level)) {
    return(list())            # hdf5r's own defaults: chunked, gzip
  }
  if (level == 0L) {
    return(list(chunk_dims = NULL))
  }
  list(gzip_level = level)
}

write_dataset <- function(parent, key, obj, scalar = FALSE, ds_args = list()) {
  dtype <- NULL
  space <- NULL
  if (is.character(obj)) {
    dtype <- H5T_STRING$new(type = "c", size = Inf)
    dtype$set_cset("UTF-8")
  }
  if (is.logical(obj)) {
    # hdf5r's default logical type is a three-value enum (FALSE/TRUE/NA), which
    # h5py reads as uint8 rather than bool. anndata rejects a non-boolean mask
    # outright, so a nullable column written with the default type cannot be
    # read back by anndata at all. Every logical reaching here is already
    # NA-free -- write_matrix() routes the rest to the nullable encoding, whose
    # own `values` and `mask` are NA-free by construction -- so the two-value
    # type loses nothing and is what other readers expect.
    dtype <- H5T_LOGICAL$new(include_NA = FALSE)
  }
  if (scalar) {
    space <- H5S$new("scalar")
  }
  do.call(
    parent$create_dataset,
    c(list(key, obj, dtype = dtype, space = space), ds_args)
  )
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

write_matrix <- function(parent, key, mat, storage_sparse_type = "csr_matrix", ds_args = list()) {
  # AnnData (>=0.8) has no nullable-string encoding, so a character vector with
  # NAs is stored as categorical, where missing values are represented by the
  # code -1. Previously NAs were coerced to the literal string "NaN", silently
  # corrupting the data.
  if (is.character(mat) && anyNA(mat)) {
    return(write_matrix(parent, key, factor(mat), storage_sparse_type, ds_args))
  }

  if (is.matrix(mat) || is.vector(mat) || is.array(mat)) {
    hasna <- anyNA(mat)
    if (hasna && is.double(mat)) {
      # FIXME: extend anndata spec to handle double NAs?
      mat[is.na(mat)] <- NaN
      hasna <- FALSE
    }

    if (!hasna) {
      dset <- write_dataset(parent, key, mat, ds_args = ds_args)
      write_attribute(dset, "encoding-type", ifelse(is.character(mat), "string-array", "array"))
      write_attribute(dset, "encoding-version", "0.2.0")
    } else {
      grp <- parent$create_group(key)
      # `values` has to be written with the NAs already filled in. Passing the
      # NA-carrying vector back to write_matrix() would land in this same branch
      # and recurse until the stack gives out -- which is what an integer column
      # holding NAs used to do (a character one is redirected to `factor` above,
      # and a double one has its NAs turned into NaN). AnnData ignores whatever
      # sits under a set mask bit, so the fill value itself is arbitrary.
      # The fill has to carry `mat`'s own type: a bare 0 is a double and would
      # silently widen an integer column to float64 on disk.
      values <- mat
      values[is.na(values)] <- as.vector(0, mode = typeof(mat))
      write_matrix(grp, "values", values, ds_args = ds_args)
      write_matrix(grp, "mask", is.na(mat), ds_args = ds_args)
      write_attribute(grp, "encoding-type", ifelse(is.logical(mat), "nullable-boolean", "nullable-integer"))
      write_attribute(grp, "encoding-version", "0.1.0")
    }
  } else if (is.factor(mat)) {
    grp <- parent$create_group(key)
    codes <- as.integer(mat)
    codes[is.na(mat)] <- 0L
    write_matrix(grp, "codes", codes - 1L, ds_args = ds_args)
    write_matrix(grp, "categories", levels(mat), ds_args = ds_args)
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
    write_dataset(grp, "indptr", mat0@p, ds_args = ds_args)
    write_dataset(grp, "data", mat0@x, ds_args = ds_args)
    write_attribute(grp, "shape", rev(dim(mat)))
    write_attribute(grp, "encoding-version", "0.1.0")
    if (is(mat0, "dgCMatrix")) {
      write_dataset(grp, "indices", mat0@i, ds_args = ds_args)
      write_attribute(grp, "encoding-type", storage_sparse_type) # "csr_matrix")
    } else {
      write_dataset(grp, "indices", mat0@j, ds_args = ds_args)
      write_attribute(grp, "encoding-type", storage_sparse_type) # "csc_matrix")
    }
  } else if (inherits(mat, "IterableMatrix")) {
    write_iterable_matrix(parent, key, mat, storage_sparse_type, ds_args)
  } else {
    stop("Writing matrices of type ", class(mat), " is not implemented: ", key)
  }
}

# Columns (cells) materialised at a time when streaming a disk-backed matrix out
# to HDF5. The point of a BPCells-backed assay is that the matrix never fits in
# memory, so the block is what bounds the peak: one block of a 30k-feature assay
# at 2k non-zeros per cell is roughly 100 MB.
.stream_block_cols <- 4096L

# Chunk length, in elements, of the extendable datasets the stream grows.
# 65536 doubles is a 512 KB chunk, large enough that gzip has something to work
# with and small enough not to inflate a short matrix.
.stream_chunk_len <- 65536L

# A dataset can only be extended if it is chunked, so the contiguous layout that
# resolve_compression() picks for compression = "none" cannot be used here.
# Keep the chunking and turn the filter off instead, which is what "none" is
# actually asking for.
resolve_stream_compression <- function(ds_args) {
  if ("chunk_dims" %in% names(ds_args) && is.null(ds_args$chunk_dims)) {
    return(list(chunk_dims = .stream_chunk_len, gzip_level = 0L))
  }
  c(ds_args, list(chunk_dims = .stream_chunk_len))
}

# An empty 1-D dataset of unlimited extent, to be grown by stream_append().
#' @import hdf5r
new_stream_dataset <- function(parent, key, dtype, ds_args) {
  do.call(
    parent$create_dataset,
    c(list(key, dtype = dtype, space = H5S$new("simple", dims = 0, maxdims = Inf)), ds_args)
  )
}

# Append to a dataset created by new_stream_dataset(). `offset` is how many
# elements it already holds; the new length is returned so callers can thread it
# through the next call. Offsets are carried as doubles rather than integers
# because the non-zero count of a matrix large enough to warrant streaming can
# exceed .Machine$integer.max.
stream_append <- function(dset, values, offset) {
  n <- length(values)
  if (n == 0L) {
    return(offset)
  }
  dset$set_extent(offset + n)
  dset[(offset + 1):(offset + n)] <- values
  offset + n
}

# Write a sparse matrix that is never held in memory in full.
#
# `dims` is the R-oriented (features x cells) shape and `block_fn(j1, j2)` returns
# columns j1..j2 of it as a dgCMatrix. Only the AnnData csr_matrix layout is
# produced: a csr_matrix is grouped by observation, a features x cells dgCMatrix
# is grouped by cell, so a block of columns appends to the output verbatim -- the
# same reinterpretation the in-memory path relies on, applied one block at a
# time. A csc_matrix would be grouped by feature, which is the opposite of the
# order a column-major source streams in.
#
# indptr is written as int64: unlike the in-memory path, which is bounded by what
# a dgCMatrix can hold, the running non-zero count here can overflow int32.
#' @import hdf5r
write_sparse_stream <- function(parent, key, dims, block_fn, ds_args = list(),
                                block_cols = .stream_block_cols) {
  grp <- parent$create_group(key)
  stream_args <- resolve_stream_compression(ds_args)

  indices <- new_stream_dataset(grp, "indices", h5types$H5T_NATIVE_INT32, stream_args)
  data <- new_stream_dataset(grp, "data", h5types$H5T_NATIVE_DOUBLE, stream_args)
  indptr <- new_stream_dataset(grp, "indptr", h5types$H5T_NATIVE_INT64, stream_args)

  n_features <- dims[1]
  n_cells <- dims[2]

  # indptr always opens with a 0 and has one entry per cell after it.
  ptr_offset <- stream_append(indptr, 0, 0)
  nnz <- 0

  if (n_cells > 0) {
    for (start in seq.int(1L, n_cells, by = block_cols)) {
      end <- min(start + block_cols - 1L, n_cells)
      block <- block_fn(start, end)
      if (!is(block, "dgCMatrix")) {
        block <- methods::as(block, "dgCMatrix")
      }
      if (nrow(block) != n_features || ncol(block) != end - start + 1L) {
        stop("Block ", start, ":", end, " of ", key, " has dimensions ",
             nrow(block), "x", ncol(block), ", expected ", n_features, "x",
             end - start + 1L, ".")
      }
      stream_append(indices, block@i, nnz)
      stream_append(data, block@x, nnz)
      # block@p is cumulative within the block and opens with a 0 that the
      # previous block already accounted for, hence [-1] and the running offset.
      ptr_offset <- stream_append(indptr, block@p[-1] + nnz, ptr_offset)
      nnz <- nnz + length(block@x)
    }
  }

  write_attribute(grp, "shape", rev(dims))
  write_attribute(grp, "encoding-type", "csr_matrix")
  write_attribute(grp, "encoding-version", "0.1.0")

  invisible(nnz)
}

# BPCells matrices are read a block of cells at a time rather than converted to
# a dgCMatrix first, which would defeat the purpose of a disk-backed assay.
# Lazily transformed matrices (what NormalizeData() leaves behind) work too: the
# transform is applied as each block is pulled.
write_iterable_matrix <- function(parent, key, mat, storage_sparse_type, ds_args) {
  if (storage_sparse_type != "csr_matrix") {
    stop("Writing a disk-backed matrix (", class(mat), ") requires ",
         "sparse.type = \"csr_matrix\": ", key, ".\n",
         "A csc_matrix is grouped by feature, the opposite of the order the ",
         "matrix streams in, so writing one would mean rewriting the whole ",
         "matrix with BPCells::transpose_storage_order() first.",
         call. = FALSE)
  }
  write_sparse_stream(
    parent, key, dim(mat),
    function(j1, j2) methods::as(mat[, j1:j2, drop = FALSE], "dgCMatrix"),
    ds_args = ds_args
  )
}

# Collect the on-disk locations the matrices of `object` read from. A lazily
# transformed BPCells matrix is a tree of operations whose leaves hold the
# paths, so the whole tree has to be walked.
matrix_backing_paths <- function(x) {
  if (!isS4(x)) {
    return(character())
  }
  paths <- character()
  for (name in methods::slotNames(class(x))) {
    value <- methods::slot(x, name)
    if (name %in% c("path", "dir") && is.character(value)) {
      paths <- c(paths, value)
    } else if (isS4(value)) {
      paths <- c(paths, matrix_backing_paths(value))
    } else if (is.list(value)) {
      paths <- c(paths, unlist(lapply(value, matrix_backing_paths), use.names = FALSE))
    }
  }
  paths
}

backing_paths <- function(object) {
  matrices <- unlist(
    lapply(object@assays, function(assay) {
      if (inherits(assay, "Assay5")) {
        return(as.list(assay@layers))
      }
      lapply(c("counts", "data", "scale.data"), function(name) {
        if (name %in% methods::slotNames(class(assay))) methods::slot(assay, name) else NULL
      })
    }),
    recursive = FALSE, use.names = FALSE
  )
  matrices <- Filter(function(m) inherits(m, "IterableMatrix"), matrices)
  if (length(matrices) == 0) {
    return(character())
  }
  paths <- unique(unlist(lapply(matrices, matrix_backing_paths), use.names = FALSE))
  normalizePath(paths, mustWork = FALSE)
}

# A disk-backed matrix keeps reading from its source for as long as the object
# is alive, and BPCells opens that source read-write whenever the file
# permissions allow it. Writing the object back over its own source would
# therefore either fail on an HDF5 lock or truncate the data mid-stream, and
# since overwrite = TRUE is the default it is easy to ask for by accident.
check_not_backing_file <- function(object, file) {
  sources <- backing_paths(object)
  if (length(sources) == 0) {
    return(invisible(NULL))
  }
  target <- normalizePath(file, mustWork = FALSE)
  if (target %in% sources) {
    stop("Cannot write to ", file, ": a disk-backed matrix in this object reads ",
         "from that same location, and writing it would destroy the data being ",
         "read.\nWrite to a different file, or bring the matrices into memory ",
         "first.", call. = FALSE)
  }
  invisible(NULL)
}

# HDF5 treats "/" as a path separator, so it cannot appear in an object name.
# Columns carrying it (e.g. "Count (cells/ul)") are renamed; names that are
# already legal are left untouched, and a rename that would collide with
# another column gets a numeric suffix.
sanitize_h5_names <- function(names) {
  dirty <- grepl("/", names, fixed = TRUE)
  if (!any(dirty)) {
    return(names)
  }
  out <- names
  for (i in which(dirty)) {
    base <- gsub("/", "_", names[i], fixed = TRUE)
    candidate <- base
    suffix <- 1L
    while (candidate %in% out[-i]) {
      candidate <- paste0(base, ".", suffix)
      suffix <- suffix + 1L
    }
    out[i] <- candidate
  }
  out
}

write_data_frame <- function(parent, key, attr_df, ds_args = list()) {
  grp <- parent$create_group(key)
  if (!is.data.frame(attr_df)) { # row names only. Creating a data.frame with duplicated row.names is not possible
    attr_df <- data.frame("_index" = attr_df, check.names = FALSE)
    attr_columns <- character()
  } else {
    attr_columns <- colnames(attr_df)
    attr_df["_index"] <- rownames(attr_df)
  }


  # The dataset name and the matching entry in "column-order" have to be derived
  # from the same sanitized name. Renaming only the dataset leaves readers that
  # iterate over column-order (anndata, mudata) unable to open the file at all.
  df_columns <- colnames(attr_df)
  h5_columns <- sanitize_h5_names(df_columns)
  renamed <- h5_columns != df_columns
  if (any(renamed)) {
    warning(
      "HDF5 does not allow '/' in names, renaming column(s): ",
      paste0(df_columns[renamed], " -> ", h5_columns[renamed], collapse = ", ")
    )
  }
  attr_columns <- h5_columns[match(attr_columns, df_columns)]

  for (i in seq_along(df_columns)) {
    col <- df_columns[i]
    h5_col <- h5_columns[i]

    # Check if the column is of (Date, POSIXct/POSIXt)
    if (inherits(attr_df[[i]], "Date") ||
      inherits(attr_df[[i]], "POSIXct") ||
      inherits(attr_df[[i]], "POSIXt")) {
      message("Column ", col, " is of datetime type.")
      attr_df[[i]] <- as.character(attr_df[[i]])
    }

    # debug: skip column with all values are NA, it's not allow
    if (all(is.na(attr_df[[i]]))) {
      attr_columns <- attr_columns[attr_columns != h5_col] # remove from column-order
      warning("Skip meta.data column: ", col, ", because of all values are NA.")
      next
    }
    write_matrix(grp, h5_col, attr_df[[i]], ds_args = ds_args)
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

# MuData records, for each modality, the 1-based index of every global
# observation/variable within that modality (0 when it is not present).
# Readers require /obsmap and /varmap to be present.
write_mod_maps <- function(h5, modalities, n_obs, var_names, ds_args = list()) {
  obsmap_group <- h5$create_group("obsmap")
  write_attribute(obsmap_group, "encoding-type", "dict")
  write_attribute(obsmap_group, "encoding-version", "0.1.0")
  varmap_group <- h5$create_group("varmap")
  write_attribute(varmap_group, "encoding-type", "dict")
  write_attribute(varmap_group, "encoding-version", "0.1.0")

  n_var <- length(unlist(var_names, use.names = FALSE))
  var_offset <- 0L
  for (mod in modalities) {
    # A Seurat object shares all of its cells across every assay.
    write_matrix(obsmap_group, mod, seq_len(n_obs), ds_args = ds_args)

    # The global var is the concatenation of the per-modality var_names,
    # so each modality occupies one contiguous block of it.
    n_mod_var <- length(var_names[[mod]])
    varmap <- integer(n_var)
    varmap[var_offset + seq_len(n_mod_var)] <- seq_len(n_mod_var)
    write_matrix(varmap_group, mod, varmap, ds_args = ds_args)
    var_offset <- var_offset + n_mod_var
  }
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
