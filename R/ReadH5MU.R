#' Read an .h5mu file and create a \code{\link{Seurat}} object.
#'
#' @param file Path to the .h5mu file.
#'
#' @return A \code{\link{Seurat}} object
#'
#' @export
ReadH5AD <- function(file) {
  # Connect to the the file
  h5 <- open_anndata(file)

  # Get metadata
  obs <- read_table(h5[["obs"]])
  var <- read_table(h5[["var"]])

  # X
  # obs and var have already been read above; passing them in keeps /obs from
  # being read once here and once again for each of obsm/obsp below.
  assay <- read_layers_to_assay(h5, obs = obs, var = var)

  # obsm
  obsm <- read_attr_m(h5, 'obs', rownames(obs))

  # varm
  varm <- read_attr_m(h5, 'var', rownames(var))

  # obsp
  obsp <- read_attr_p(h5, 'obs', rownames(obs))

  # If there are var pairs, there's no place to store it
  # in the Seurat object
  # var_pairs <- read_attr_p(h5, 'var')
  var_pairs_names <- c()
  if ("varp" %in% names(h5))
    var_pairs_names <- names(h5[["varp"]])
  if (!is.null(var_pairs_names) && !isFALSE(var_pairs_names) && length(var_pairs_names) > 0)
    missing_on_read("/varp", "pairwise annotation of variables")

  # Create a Seurat object
  # If read from .h5mu modality, give an assay name
  path_fragments <- strsplit(file, "\\.h5mu")[[1]]
  assay_name <- "RNA"  # Seurat's own default
  if (length(path_fragments) == 2) {
    mod_path_fragments <- strsplit(path_fragments[2], "\\/")[[1]]
    assay_name <- mod_path_fragments[length(mod_path_fragments)]
  }
  if (skip_calcn(assay, assay_name, colnames(obs))) {
    srt <- without_calcn(Seurat::CreateSeuratObject(assay, assay = assay_name))
  } else {
    srt <- Seurat::CreateSeuratObject(assay, assay = assay_name)
  }

  # Specify highly variable features
  if ("highly_variable" %in% colnames(var) && is.logical(var$highly_variable)) {
    Seurat::VariableFeatures(srt[[assay_name]]) <-
      rownames(srt[[assay_name]])[var$highly_variable]
  }

  # Add metadata
  srt@meta.data <- add_meta_data(srt@meta.data, obs)

  # NOTE: feature metadata is attached by read_layers_to_assay(). It used to be
  # cbind()ed on again here, which duplicated every column of /var.

  # Add embeddings
  for (emb in names(obsm)) {
    emb_name <- gsub('X_', '', emb)

    maybe_loadings <- matrix()
    if (emb %in% names(OBSM2VARM)) {
      varm_key = OBSM2VARM[[emb]]
      if (varm_key %in% names(varm)) {
        maybe_loadings <- varm[[varm_key]]
      }
    }

    emb_stdev <- numeric()
    if ("uns" %in% names(h5)) {
      if (emb_name %in% names(h5[["uns"]])) {
        if ("variance" %in% names(h5[["uns"]][[emb_name]])) {
          emb_stdev <- sqrt(h5[["uns"]][[emb_name]][["variance"]]$read())
        }
      }
    }

    srt[[emb_name]] <- Seurat::CreateDimReducObject(
      embeddings = obsm[[emb]][rownames(obs),,drop=FALSE],
      loadings = maybe_loadings,
      key = paste0(emb_name, "_"),
      assay = Seurat::DefaultAssay(srt),
      stdev = emb_stdev
    )
  }

  # Add graphs
  srt@graphs <- lapply(obsp, Seurat::as.Graph)

  # Close the connection
  h5$close()

  srt
}

#' Create a \code{Seurat} object from .h5mu file contents
#'
#' @param file Path to the .h5mu file
#'
#' @import hdf5r Matrix Seurat
#' @importFrom utils hasName
#'
#' @return A \code{Seurat} object
#' '
#' @export ReadH5MU
ReadH5MU <- function(file) {
  # Connect to the the file
  h5 <- open_and_check_mudata(file)

  # Get assays (modalities)
  assays <- h5[["mod"]]$names
  if ("mod-order" %in% names(h5attributes(h5[["mod"]]))) {
    modorder <- h5attributes(h5[["mod"]])$`mod-order`
    if (all(assays %in% modorder)) {
      modorder <- modorder[modorder %in% assays]
      if (!any(duplicated(modorder))) {
        assays <- modorder
      }
    }
  }

  # Get global metadata
  metadata <- read_table(h5[["obs"]])

  # NOTE: there's no global feature metadata in the Seurat object
  ft_metadata <- tryCatch({
      read_table(h5[["var"]])
    },
    error = function(err) {
      warning(err)
      read_table(h5[["var"]], set_index = FALSE)
    }
  )

  if (ncol(ft_metadata) > 0)
    missing_on_read("/var", paste0("global variables metadata (", paste(colnames(ft_metadata), collapse = ", "), ")"))

  # Get (multimodal) embeddings
  embeddings <- read_attr_m(h5, 'obs', rownames(metadata))
  # If obs->mod mappings are in the file, dismiss them
  embeddings <- embeddings[!names(embeddings) %in% assays]

  # Get (multimodal) loadings
  # NOTE: features can be set as row names only if they are unique
  loadings <- read_attr_m(h5, 'var', rownames(ft_metadata))

  # Get obs pairs
  obs_pairs <- read_attr_p(h5, 'obs', rownames(metadata))

  # If there are var pairs, there's no place to store it
  # in the Seurat object
  var_pairs_names <- c()
  if ("varp" %in% names(h5))
    var_pairs_names <- names(h5[["varp"]])
  if (!is.null(var_pairs_names) && !isFALSE(var_pairs_names) && length(var_pairs_names) > 0)
    missing_on_read("/varp", "pairwise annotation of variables")

  # mod/.../obs and mod/.../var
  # These are read up front because the assay, obsm, varm and obsp of a modality
  # all need them; reading them once here rather than inside each of those steps
  # avoids re-reading every modality's /obs three more times.
  mod_obs <- lapply(assays, function(mod) {
    read_table(h5[['mod']][[mod]][['obs']])
  })
  names(mod_obs) <- assays

  mod_var <- lapply(assays, function(mod) {
    read_table(h5[['mod']][[mod]][['var']])
  })
  names(mod_var) <- assays

  # mod/.../X, raw, and layers
  modalities <- lapply(assays, function(mod) {
   read_layers_to_assay(h5[['mod']][[mod]], mod, obs = mod_obs[[mod]], var = mod_var[[mod]])
  })
  names(modalities) <- assays

  # mod/.../obsm
  mod_obsm <- lapply(assays, function(mod) {
    read_attr_m(h5[['mod']][[mod]], 'obs', rownames(mod_obs[[mod]]))
  })
  names(mod_obsm) <- assays

  # mod/.../varm
  mod_varm <- lapply(assays, function(mod) {
    read_attr_m(h5[['mod']][[mod]], 'var', rownames(mod_var[[mod]]))
  })
  names(mod_varm) <- assays

  # mod/.../obsp
  mod_obsp <- lapply(assays, function(mod) {
    read_attr_p(h5[['mod']][[mod]], 'obs', rownames(mod_obs[[mod]]))
  })
  names(mod_obsp) <- assays

  # If there are var pairs in individual modalities,
  # there's no place to store it in the Seurat object.
  for (mod in assays) {
    if ("varp" %in% names(h5[['mod']][[mod]])) {
      if (length(h5[['mod']][[mod]][['varp']]) > 0) {
        missing_on_read(paste0("/mod", mod, "/varp"), "pairwise annotation of variables")
      }
    }
  }

  var_pairs_names <- c()
  if ("varp" %in% names(h5))
    var_pairs_names <- names(h5[["varp"]])
  if (!is.null(var_pairs_names) && !isFALSE(var_pairs_names) && length(var_pairs_names) > 0)
    missing_on_read("/varp", "pairwise annotation of variables")

  # Only common observations can be read
  obs_names <- Reduce(intersect, lapply(modalities, colnames))
  mods_n_obs <- unique(vapply(modalities, ncol, 1))
  if (length(mods_n_obs) > 1 || mods_n_obs[1] != length(obs_names)) {
    warning("Only the intersection of observations (samples) is loaded. Observations that are not present in all the modalities (assays) are discarded.")
  }

  # Create a Seurat object
  # Only CreateSeuratObject() recomputes nCount/nFeature, and it does so for the
  # first modality only; assigning the remaining assays below does not.
  first_assay <- subset_cells(modalities[[1]], obs_names)
  merged_meta_columns <- c(colnames(metadata), unlist(lapply(mod_obs, colnames), use.names = FALSE))
  if (skip_calcn(first_assay, names(modalities)[1], merged_meta_columns)) {
    srt <- without_calcn(Seurat::CreateSeuratObject(first_assay, assay = names(modalities)[1]))
  } else {
    srt <- Seurat::CreateSeuratObject(first_assay, assay = names(modalities)[1])
  }
  # NOTE: [-1], not [2:length()], which yields c(NA, ..) for one modality
  for (modality in names(modalities)[-1]) {
    srt[[modality]] <- subset_cells(modalities[[modality]], obs_names)
  }

  # Global /obs is where a Seurat object's meta.data is written, so it has to be
  # attached; previously it was read and then discarded, losing all of it.
  if (ncol(metadata) > 0 && all(obs_names %in% rownames(metadata))) {
    srt@meta.data <- add_meta_data(srt@meta.data, metadata[obs_names, , drop = FALSE])
  } else if (ncol(metadata) > 0) {
    warning("Global /obs could not be matched to the observation names and was not loaded.")
  }

  # Metadata, features metadata, and variable features
  for (modality in names(modalities)) {
    # Append modality metadata
    srt@meta.data <- add_meta_data(srt@meta.data, mod_obs[[modality]][obs_names, , drop = FALSE])

    metafeatures <- srt[[modality]][[]]

    # Specify highly variable features
    if ("highly_variable" %in% colnames(metafeatures) &&
        is.logical(metafeatures$highly_variable)) {
      Seurat::VariableFeatures(srt[[modality]]) <-
        rownames(metafeatures)[metafeatures$highly_variable]
    }
  }

  # Add joint embeddings
  for (emb in names(embeddings)) {
    emb_name <- toupper(gsub('X_', '', emb))

    maybe_loadings <- matrix()
    varm_key <- emb_name
    if (emb %in% names(OBSM2VARM)) {
      varm_key <- OBSM2VARM[[emb]]
    }
    if (hasName(loadings, varm_key)) {
      maybe_loadings <- loadings[[varm_key]]
    }

    emb_stdev <- numeric()
    if ("uns" %in% names(h5)) {
      if (emb_name %in% names(h5[["uns"]])) {
        if ("variance" %in% names(h5[["uns"]][[emb_name]])) {
          emb_stdev <- sqrt(h5[["uns"]][[emb_name]][["variance"]]$read())
        }
      }
    }

    srt[[emb_name]] <- Seurat::CreateDimReducObject(
      embeddings = embeddings[[emb]][obs_names,,drop=FALSE],
      loadings = maybe_loadings,
      key = paste0(emb_name, "_"),
      stdev = emb_stdev,
      assay = Seurat::DefaultAssay(srt),  # this is not true but an existing assay must be provided
    )
  }

  # Do embeddings across modalities have unique names?
  # If each modality has e.g. X_pca, we will need to harmonise the names
  # by prepending modality name: e.g. RNAPCA.
  unique_emb <- FALSE
  all_embeddings <- c(names(embeddings), unlist(lapply(mod_obsm, function(obsm) names(obsm))))
  if (!any(duplicated(all_embeddings))) {
    unique_emb <- TRUE
  }

  # Add modality-specific embeddings
  for (mod in names(mod_obsm)) {
    mod_embeddings <- mod_obsm[[mod]]
    for (emb in names(mod_embeddings)) {
      emb_name <- gsub('X_', '', emb)
      if (!startsWith(emb_name, mod))
        emb_name <- toupper(emb_name)
      modemb_name <- emb_name
      if (!unique_emb)
        modemb_name <- paste(mod, emb_name, sep = "")
      # Embeddings keys will have to follow the format alphanumericcharacters_, e.g. RNAPCA_.

      maybe_loadings <- matrix()
      varm_key <- emb_name
      if (emb %in% names(OBSM2VARM))
        varm_key = OBSM2VARM[[emb]]
      if (hasName(mod_varm[[mod]], varm_key))
        maybe_loadings <- mod_varm[[mod]][[varm_key]]

      emb_stdev <- numeric()
      h5_mod <- h5[["mod"]][[mod]]
      if ("uns" %in% names(h5_mod)) {
        if (emb_name %in% names(h5_mod[["uns"]])) {
          if ("variance" %in% names(h5_mod[["uns"]][[emb_name]])) {
            emb_stdev <- sqrt(h5_mod[["uns"]][[emb_name]][["variance"]]$read())
          }
        }
      }

      srt[[modemb_name]] <- Seurat::CreateDimReducObject(
        embeddings = mod_embeddings[[emb]][obs_names,,drop=FALSE],
        loadings = maybe_loadings,
        key = paste0(modemb_name, "_"),
        stdev = emb_stdev,
        assay = mod,
      )
    }
  }

  # Add graphs

  # Only take into account common observations
  if (length(obs_pairs) > 0) {
    srt@graphs <- lapply(obs_pairs, function(graph) {
      graph[obs_names,obs_names,drop=FALSE]
    })
    names(srt@graphs) <- names(obs_pairs)
  }

  for (mod in names(mod_obsp)) {
    for (graph in names(mod_obsp[[mod]])) {
      graph_name <- graph
      # /mod/RNA/obsp/distances -> @graphs$RNA_distances
      if (graph %in% names(srt@graphs)) {
        graph_name <- paste(mod, graph, sep = "_")
      }
      srt@graphs[[graph_name]] <- Seurat::as.Graph(mod_obsp[[mod]][[graph]][obs_names, obs_names])
      srt@graphs[[graph_name]]@assay.used <- mod
    }
  }


  # Close the connection
  h5$close_all()

  srt
}
