#' @title MovidaData
#' @description An R6 class for managing and analyzing multi-omics data.
#' @details The `MovidaModel` class provides methods for initializing, processing,
#' and retrieving relationships and data from proteomics, transcriptomics,
#' and metabolomics experiments. It supports parallel processing and
#' includes functionality for filtering and retrieving data based on
#' specific criteria.
#'
#' @section Imports:
#' This class requires the following functions and packages:
#' - `R6::R6Class`
#' - `SummarizedExperiment::SummarizedExperiment`, `rowData`, `colData`, `assay`
#' - `future::plan`, `future::multisession`, `future::sequential`
#' - `parallel::detectCores`
#' - `readRDS`, `saveRDS`
#' - Internal functions: `build_inchi_relationships`, `build_uniprot_to_ensembl`, `check_movida_list`, `check_ensembl`, `check_inchi`, `check_uniprot`
#'
#' @importFrom R6 R6Class
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData assay
#' @importFrom future plan multisession sequential
#' @importFrom parallel detectCores
#' @importFrom methods is
#' @section Private Fields:
#' \describe{
#'   \item{dde_prot}{SummarizedExperiment object for proteomics data.}
#'   \item{dde_trans}{SummarizedExperiment object for transcriptomics data.}
#'   \item{dde_metabo}{SummarizedExperiment object for metabolomics data.}
#'   \item{inchi_to_ensembl}{Relationships between InChI and Ensembl IDs.}
#'   \item{inchi_to_uniprot}{Relationships between InChI and UniProt IDs.}
#'   \item{uniprot_to_ensembl}{Relationships between UniProt and Ensembl IDs.}
#'   \item{organism}{Organism information for the data.}
#' }
#'
#' @section Public Methods:
#' The class provides a suite of methods for interacting with the data model.
#'
#' - `initialize(movida_list)`: Initializes the object with input data for proteomics, transcriptomics,
#'   and metabolomics, along with associated differential analysis and enrichment results.
#'
#' - `build_relationships()`: Constructs feature mappings between omics layers using InChI, UniProt,
#'   and Ensembl identifiers, with support for parallel processing.
#'
#' - `load_relationships(folder_path)`: Loads previously saved relationship data from disk.
#'
#' - `save_relationships(folder_path)`: Saves current feature relationship data to disk for reuse.
#'
#' - `get_contrasts()`: Returns a list of available contrasts (comparisons) in the analysis.
#'
#' - `getDEA(target, contrast, FDRpvalue = NULL, FDRadj = NULL)`: Retrieves differential expression results
#'   for the given data target and contrast, with optional FDR-based filtering.
#'
#' - `getFEA(target, contrast, FDRpvalue = NULL, FDRadj = NULL)`: Retrieves functional enrichment analysis results
#'   for the specified target and contrast, with optional filtering.
#'
#' - `get_values_all(target, return_se = FALSE)`: Returns all assay data (or `SummarizedExperiment` object)
#'   for the selected data target.
#'
#' - `get_valuesfeatures(features, source = NULL)`: Retrieves expression data for specific features.
#'
#' - `get_values_group(groups, source)`: Returns data for specified experimental groups.
#'
#' - `get_values_samples(samples, source)`: Returns data for specific sample IDs.
#'
#' - `get_related_features(feature, target)`: Maps a feature identifier (e.g., InChI, UniProt, Ensembl)
#'   to its related identifiers in other omics layers.
#'
#' - `get_inchi_to_ensembl()`, `get_inchi_to_uniprot()`, `get_uniprot_to_ensembl()`: Accessors for relationship mappings
#'   generated during the `build_relationships()` phase.
#'
#'
#' @examples
#' # Example usage:
#' movida_list <- list(
#'   results_prot = list(),
#'   results_trans = list(),
#'   results_metabo = list(),
#'   dde_prot = NULL,
#'   dde_trans = NULL,
#'   dde_metabo = NULL,
#'   organism = "Hs"
#' )
#' model <- MovidaModel$new(movida_list)
#' @docType class
#' @export
MovidaModel <- R6Class("MovidaModel",
  private = list(
    dde_prot = NULL,
    dde_trans = NULL,
    dde_metabo = NULL,
    inchi_to_ensembl = NULL,
    inchi_to_uniprot = NULL,
    uniprot_to_ensembl = NULL,
    organism = NULL,
    get_dde_object = function(source) {
      # Validate that 'source' is a single string
      if (!is.character(source) || length(source) != 1) {
        warning("get_dde_object: Argument 'source' must be a single string.")
        return(NULL)
      }

      # Select the appropriate SummarizedExperiment object based on the source
      if (source == "proteomics") {
        dde <- private$dde_prot
      } else if (source == "transcriptomics") {
        dde <- private$dde_trans
      } else if (source == "metabolomics") {
        dde <- private$dde_metabo
      } else {
        stop("get_dde_object: Invalid source type.")
      }
      return(dde)
    }
  ),
  public = list(
    #' @description Initialize the MovidaModel object with data.
    #' @param movida_list A list containing the data for proteomics, transcriptomics, and metabolomics.
    initialize = function(movida_list) {
      check_movida_list_dde(movida_list)

      # Initialize private variables with data from the provided movida_list
      private$dde_prot <- movida_list$dde_prot
      private$dde_trans <- movida_list$dde_trans
      private$dde_metabo <- movida_list$dde_metabo
      private$organism <- movida_list$organism

      if (!is.null(movida_list$metadata)) {
        if (!check_metadata(movida_list$metadata)) {
          warning("Metadata does not match the data. Please check your metadata.")
        } else {
          private$metadata <- movida_list$metadata

          colData(private$dde_prot) <- S4Arrays::DataFrame(movida_list$metadata)
          colData(private$dde_trans) <- S4Arrays::DataFrame(movida_list$metadata)
          colData(private$dde_metabo) <- S4Arrays::DataFrame(movida_list$metadata)
        }
      }
    },
    #' @description Get the SummarizedExperiment object for a specific source.
    #'
    #' @param source A character string indicating the source type: "proteomics", "transcriptomics", or "metabolomics".
    #'
    #' @return Returns the corresponding SummarizedExperiment object.
    get_dde_object_exposed = function(source) {
      private$get_dde_object(source)
    },
    #' @description Funtion to build the relationships datasetss
    build_relationships = function() {
      # Use parallel processing for building relationships!! used for inchi relationships
      if (inherits(future::plan(), "multisession")) {
        future::plan(future::sequential) # Stop the existing multisession plan
        warning("Parallel processing was already set to multisession. Restarting it.")
      }
      future::plan(future::multisession, workers = parallel::detectCores(logical = FALSE)) # Set the number of workers as needed

      # Build inchi relationships if metabolomics data is available
      if (!is.null(private$dde_metabo)) {
        # Build inchi to Ensembl relationships
        message("Creating relationships dataset inchi_to_uniprot...")
        if (!is.null(private$dde_prot)) {
          private$inchi_to_ensembl <- build_inchi_relationships(rowData(private$dde_metabo), rowData(private$dde_prot), get_ensembl = FALSE, get_uniprot = TRUE)
        }
        message("Done!")

        # Build inchi to UniProt relationships
        message("Creating relationships dataset inchi_to_ensembl...")
        if (!is.null(private$dde_trans)) {
          private$inchi_to_uniprot <- build_inchi_relationships(rowData(private$dde_metabo), rowData(private$dde_trans), get_ensembl = TRUE, get_uniprot = FALSE)
        }
        message("Done!")
      }

      # Build UniProt-to-Ensembl relationships if both proteomics and transcriptomics data are available
      if (!is.null(private$dde_prot) && !is.null(private$dde_trans)) {
        private$uniprot_to_ensembl <- build_uniprot_to_ensembl(rowData(private$dde_prot), rowData(private$dde_trans), private$organism)
      }
    },

    #' @description Unify all the contrasts from the movida_list components.
    #'
    #' @return A character vector of unique contrast names found across the specified elements of `movida_list`.
    get_contrasts = function() {
      return(
        unique(c(
          suppressWarnings(DEANames_wrapper(private$dde_metabo)),
          suppressWarnings(DEANames_wrapper(private$dde_prot)),
          suppressWarnings(DEANames_wrapper(private$dde_trans))
        ))
      )
    },

    #' @description Load Relationships from Folder
    #'
    #' @param folder_path Character string specifying the path to the folder containing relationship data files.
    #'
    #' @return Returns a data.frame containing the loaded relationships.
    #'
    #' @examples
    #' relationships <- load_relationships("path/to/relationships_folder")
    load_relationships = function(folder_path) {
      if (!is.character(folder_path) || length(folder_path) != 1) {
        warning("Argument 'folder_path' must be a single string.")
        return(NULL)
      }
      if (!dir.exists(folder_path)) {
        warning("The specified folder does not exist.")
        return(NULL)
      }

      private$inchi_to_ensembl <- readRDS(file.path(folder_path, "inchi_to_ensembl.rds"))
      private$inchi_to_uniprot <- readRDS(file.path(folder_path, "inchi_to_uniprot.rds"))
      private$uniprot_to_ensembl <- readRDS(file.path(folder_path, "uniprot_to_ensembl.rds"))
      return(invisible(TRUE))
    },

    #' @description Save relationships between features to disk.
    #'
    #' @param folder_path Path to the folder where relationships will be saved.
    #'
    #' @return Returns nothing. Saves relationship data as RDS files.
    save_relationships = function(folder_path) {
      if (!is.character(folder_path) || length(folder_path) != 1) {
        warning("Argument 'folder_path' must be a single string.")
      }
      if (!dir.exists(folder_path)) {
        dir.create(folder_path, recursive = TRUE)
      }

      saveRDS(private$inchi_to_ensembl, file.path(folder_path, "inchi_to_ensembl.rds"))
      saveRDS(private$inchi_to_uniprot, file.path(folder_path, "inchi_to_uniprot.rds"))
      saveRDS(private$uniprot_to_ensembl, file.path(folder_path, "uniprot_to_ensembl.rds"))
    },

    #' @description Get mapping from InChI to Ensembl.
    #'
    #' @return Returns a data frame with InChI to Ensembl relationships.
    get_inchi_to_ensembl = function() {
      return(private$inchi_to_ensembl)
    },

    #' @description Get mapping from InChI to UniProt.
    #'
    #' @return Returns a data frame with InChI to UniProt relationships.
    get_inchi_to_uniprot = function() {
      return(private$inchi_to_uniprot)
    },

    #' @description Get mapping from UniProt to Ensembl.
    #'
    #' @return Returns a data frame with UniProt to Ensembl relationships.
    get_uniprot_to_ensembl = function() {
      return(private$uniprot_to_ensembl)
    },
    # #' @description Get differential expression analysis results.
    # #'
    # #' @param source Data source: "proteomics", "transcriptomics", or "metabolomics".
    # #' @param contrast Contrast name.
    # #' @param FDRpvalue Optional FDR p-value threshold.
    # #' @param FDRadj Optional FDR adjusted value threshold.
    # #'
    # #' @return Returns a data frame of DEA results filtered by FDR thresholds.
    # get_objects = function() {
    #   return(names(private$object_list))
    # },

    #' @description Get differential expression analysis results.
    #'
    #' @param source Data source: "proteomics", "transcriptomics", or "metabolomics".
    #' @param contrast Contrast name.
    #' @param FDRpvalue Optional FDR p-value threshold.
    #' @param FDRadj Optional FDR adjusted value threshold.
    #'
    #' @return Returns a data frame of DEA results filtered by FDR thresholds.
    getDEA = function(source, contrast, FDRpvalue = NULL, FDRadj = NULL) {
      if (!is.character(contrast) || length(contrast) != 1) {
        warning("Argument 'contrast' must be a single string.")
        return(NULL)
      }

      data <- DEA_wrapper(private$get_dde_object(source), contrast, format = "original")

      if (is.null(data)) {
        warning("No data found for the specified source and contrast.")
        return(NULL)
      }
      data <- as.data.frame(data)
      if (!is.null(FDRpvalue)) {
        data <- data[data$FDRpvalue <= FDRpvalue, ]
      }
      if (!is.null(FDRadj)) {
        data <- data[data$FDRadj <= FDRadj, ]
        if (!is.null(FDRpvalue)) {
          warning("Both FDRpvalue and FDRadj are set. Filtering will be applied based on both criteria.")
        }
      }
      return(data)
    },

    #' @description Get functional enrichment analysis results.
    #'
    #' @param source Data source: "proteomics", "transcriptomics", or "metabolomics".
    #' @param contrast Contrast name.
    #' @param FDRpvalue Optional FDR p-value threshold.
    #' @param FDRadj Optional FDR adjusted value threshold.
    #'
    #' @return Returns a data frame of FEA results filtered by FDR thresholds.
    getFEA = function(source, contrast, FDRpvalue = NULL, FDRadj = NULL) {
      # if (!is.character(source) || length(source) != 1) {
      #   warning("Argument 'source' must be a single string.")
      #   return(NULL)
      # }

      if (!is.character(contrast) || length(contrast) != 1) {
        warning("Argument 'contrast' must be a single string.")
        return(NULL)
      }
      data <- FEA_wrapper(private$get_dde_object(source), contrast, format = "original")

      if (is.null(data)) {
        warning("No data found for the specified source and contrast.")
        return(NULL)
      }

      data <- as.data.frame(data)
      if (!is.null(FDRpvalue)) {
        data <- data[data$FDRpvalue <= FDRpvalue, ]
      }
      if (!is.null(FDRadj)) {
        data <- data[data$FDRadj <= FDRadj, ]
        if (!is.null(FDRpvalue)) {
          warning("Both FDRpvalue and FDRadj are set. Filtering will be applied based on FDRadj")
        }
      }
      return(data)
    },

    #' @description Get all values from the specified source.
    #'
    #' @param source Data source: "proteomics", "transcriptomics", or "metabolomics".
    #' @param return_se If TRUE, returns the SummarizedExperiment object; otherwise returns assay data.
    #'
    #' @return Returns assay data or SummarizedExperiment object.
    get_values_all = function(source, return_se = FALSE) {
      # Validate that 'source' is a single string
      # if (!is.character(source) || length(source) != 1) {
      #   warning("Argument 'source' must be a single string.")
      #   return(NULL)
      # }

      # Select the appropriate SummarizedExperiment object based on the source
      dde <- private$get_dde_object(source)

      # Return the assay data or the SummarizedExperiment object based on 'return_se'
      if (return_se) {
        return(dde) # Return the assay data (matrix of values)
      } else {
        return(assay(dde)) # Return the SummarizedExperiment object
      }
    },

    #' @description Get related features across omics types.
    #'
    #' @param feature feature identifier (Ensembl, InChI, or UniProt).
    #' @param target Target omics type: "proteomics", "transcriptomics", or "metabolomics".
    #'
    #' @return Returns related features in the target omics type.
    get_related_features = function(feature, source, target) {
      if (!is.character(c(feature))) {
        message("Argument 'features' must be a vector of strings.")
        return(NULL)
      }

      # Validate that the features are of a recognized type (Ensembl, inchi, or UniProt)
      if (!suppressWarnings(check_is_valid_feature(c(feature)))) {
        message("Error:get_related_features, features must be of type Ensembl, inchi, or UniProt.")
        return(NULL)
      }

      # Check if 'target' is a single string
      if (!is.character(target) || length(target) != 1) {
        message("Argument 'target' must be a single string.")
        return(NULL)
      }

      # Return the related features based on the target type
      # Allow 'source' and 'target' to be matched by partial string (e.g., if source contains "transcript", etc.)
      src <- tolower(source)
      tgt <- tolower(target)

      if (grepl("transcript", src)) {
        if (grepl("proteom", tgt)) {
          return(private$uniprot_to_ensembl[private$uniprot_to_ensembl$ENSEMBL %in% feature, , drop = FALSE]$UNIPROT)
        } else if (grepl("transcript", tgt)) {
          return(c(feature)) # Ensembl to Ensembl is a direct match
        } else if (grepl("metabolom", tgt)) {
          return(private$inchi_to_ensembl[private$inchi_to_ensembl$ENSEMBL %in% feature, , drop = FALSE]$INCHIKEY)
        }
      } else if (grepl("metabolom", src)) {
        if (grepl("proteom", tgt)) {
          return(private$inchi_to_uniprot[private$inchi_to_uniprot$INCHIKEY %in% feature, , drop = FALSE]$UNIPROT)
        } else if (grepl("transcript", tgt)) {
          return(private$inchi_to_ensembl[private$inchi_to_ensembl$INCHIKEY %in% feature, , drop = FALSE]$ENSEMBL)
        } else if (grepl("metabolom", tgt)) {
          return(c(feature)) # inchi to inchi is a direct match
        }
      } else if (grepl("proteom", src)) {
        if (grepl("proteom", tgt)) {
          return(c(feature)) # UniProt to UniProt is a direct match
        } else if (grepl("transcript", tgt)) {
          return(private$uniprot_to_ensembl[private$uniprot_to_ensembl$UNIPROT %in% feature, , drop = FALSE]$ENSEMBL)
        } else if (grepl("metabolom", tgt)) {
          return(private$inchi_to_uniprot[private$inchi_to_uniprot$UNIPROT %in% feature, , drop = FALSE]$INCHIKEY)
        }
      }
    },

    #' @description Get values for specified features.
    #'
    #' @param features Vector of feature identifiers.
    #' @param source Data source "proteomics", "transcriptomics", or "metabolomics".
    #' @param return_se If TRUE, returns the SummarizedExperiment object, otherwise returns assay data.
    #'
    #' @return Returns assay data for the specified features.
    get_values_features = function(features, source, return_se = FALSE) {
      # Check if 'features' is a vector of strings
      if (!is.character(features) || !is.vector(features)) {
        warning("Argument 'features' must be a vector of strings.")
        return(NULL)
      }

      # Select the appropriate SummarizedExperiment object based on the source
      dde_source <- private$get_dde_object(source)

      # Handle errors
      if (!all(features %in% rownames(dde_source))) {
        if (length(features) == 0) {
          warning("No features provided. Returning NULL.")
          return(NULL)
        }
        subset <- features[features %in% rownames(dde_source)]
        if (length(subset) == 0) {
          warning("No features found in the dataset.")
          return(NULL)
        }
        if (length(subset) < length(features)) {
          warning("Some features are not present in the dataset. Returning only available features.")
          features <- subset
        }
      }
      if (return_se) {
        return(dde_source[features, ]) # Return the SummarizedExperiment object
      } else {
        return(assay(dde_source)[features, , drop = FALSE])
      }
    },

    #' @description Get features associated with a pathway.
    #'
    #' @param pathway Pathway name.
    #' @param contrast Contrast name.
    #' @param source Data source "proteomics", "transcriptomics", or "metabolomics".
    #'
    #' @return Returns a vector of feature identifiers in the pathway.
    get_pathway_features = function(pathway, contrast, source) {
      # Validate that 'pathway' is a single string
      if (!is.character(pathway) || length(pathway) != 1 || suppressWarnings(check_goterm(pathway))) {
        warning("Argument 'pathway' must be a single string.")
        return(NULL)
      }

      enrich_res <- self$getFEA(source, contrast)

      if (is.null(enrich_res)) {
        warning("No enrichment results found for the specified pathway and contrast.")
        return(NULL)
      }

      # Check that it respects DeeDee standard!!
      # thisset_name <- enrich_res[pathway, "Term"] # substitute to have the standard in deedee
      thisset_members <- unlist(strsplit(enrich_res[pathway, "genes"], ","))
      # thisset_members_ens <- rowData(se)$ENSEMBL[match(thisset_members, rowData(se)$SYMBOL)]

      # if (is.na(thisset_name) || is.na(thisset_members_ens)) {
      #   warning("Pathway not found in enrichment results or no members associated with it. Check logic.")
      #   return(c(""))
      # }

      return(thisset_members)
    },

    #' @description Get values for a subset of samples based on metadata.
    #'
    #' @param subset Vector of group names to subset.
    #' @param source Data source: "proteomics", "transcriptomics", or "metabolomics".
    #' @param column Metadata column to use for subsetting.
    #' @param return_se If TRUE, returns the SummarizedExperiment object; otherwise returns assay data.
    #'
    #' @return Returns assay data or SummarizedExperiment object for the subset.
    get_values_subset_metadata = function(subset, source, column = "group", return_se = FALSE) {
      # Check if 'groups' is a vector of strings

      if (is.null(subset)) {
        return(NULL)
      }

      if (!is.character(subset) || !is.vector(subset) || length(subset) == 0) {
        warning("Argument 'subset' must be a vector of strings.")
        return(NULL)
      }

      if (!is.character(column) || length(column) != 1) {
        warning("Argument 'column' must be a single string.")
        return(NULL)
      }

      # Select the appropriate SummarizedExperiment object based on the source
      dde_source <- private$get_dde_object(source)

      # Check that the column exists in colData(dde_source)
      if (!column %in% colnames(colData(dde_source))) {
        warning(paste("Column", column, "not found in metadata."))
        return(NULL)
      }

      # Check that all groups in 'groups' are present in colData(dde_source)$group
      if (!all(subset %in% colData(dde_source)[[column]])) {
        subset_subset <- subset[subset %in% colData(dde_source)[[column]]]
        if (length(subset_subset) == 0) {
          warning("No subset found in the metadata. Returning NULL.")
          return(NULL)
        }
        if (length(subset_subset) < length(subset)) {
          warning("Some features are not present in the dataset. Returning only available features.")
          subset <- subset_subset
        }
      }

      # Subset dde_source based on the specified groups
      subset_indices <- colData(dde_source)[[column]] %in% subset
      if (return_se) {
        return(dde_source[, subset_indices]) # Return the SummarizedExperiment object
      } else {
        return(assay(dde_source)[, subset_indices, drop = FALSE]) # Return the assay data (matrix of values)
      }
    },

    #' @description Get values for specified samples.
    #'
    #' @param samples Vector of sample names.
    #' @param source String: "proteomics", "transcriptomics", or "metabolomics".
    #' @param return_se Bool: whether to return a SummarizedExperiment object or just the assay data.
    #'
    #' @return Returns SummarizedExperiment object for the specified samples.
    get_values_samples = function(samples, source, return_se = FALSE) {
      # Check if 'samples' is a vector of strings
      if (!is.character(samples) || !is.vector(samples)) {
        warning("Argument 'samples' must be a vector of strings.")
        return(NULL)
      }

      # Select the appropriate SummarizedExperiment object based on the source
      dde_source <- private$get_dde_object(source)

      # Check that all samples in 'samples' are present in colnames(dde_source)
      if (!all(samples %in% colnames(dde_source))) {
        subset <- samples[samples %in% colnames(dde_source)]
        if (length(subset) == 0) {
          warning("No subset found in the metadata. Returning NULL.")
          return(NULL)
        }
        if (length(subset) < length(samples)) {
          warning("Some features are not present in the dataset. Returning only available features.")
          samples <- subset
        }
      }


      # Subset dde_source based on the specified samples
      subset_indices <- colnames(dde_source) %in% samples
      if (return_se) {
        return(dde_source[, subset_indices]) # Return the SummarizedExperiment object
      } else {
        return(assay(dde_source)[, subset_indices, drop = FALSE]) # Return the assay data (matrix of values)
      }
    },

    #' @description Get list of features for a source.
    #'
    #' @param source Data source: "proteomics", "transcriptomics", or "metabolomics".
    #'
    #' @return Returns a data frame of feature identifiers.
    get_features_all = function(source) {
      return(rownames(private$get_dde_object(source)))
    },

    #' @description Get metadata column names for a source.
    #'
    #' @param source Data source: "proteomics", "transcriptomics", or "metabolomics".
    #'
    #' @return Returns a vector of metadata column names.
    get_metadata_columns = function(source) {

      colnames <- colnames(private$get_dde_object(source))
      colData_obj <- colData(private$get_dde_object(source))

      if (is.null(colnames)) {
        message("Get metadata columns: No metadata columns found for the specified source.")
        return(NULL)
      }


      valid_cols <- !apply(as.data.frame(colData_obj), 2, function(col) all(is.na(col) | is.null(col)))
      colnames <- colnames(colData_obj)[valid_cols]

      return(colnames)
    },

    #' @description Get metadata for a source.
    #'
    #' @param source Data source: "proteomics", "transcriptomics", or "metabolomics".
    #'
    #' @return Returns metadata for the specified source.
    get_metadata = function(source) {
      return(data.frame(colData(private$get_dde_object(source))))
    },
    #' @description Get the possible sources
    #'
    #' @param source Data source: "proteomics", "transcriptomics", or "metabolomics".
    #'
    #' @return Returns the sources available in the model.
    get_sources = function() {
      sources <- c()
      if (!is.null(private$dde_trans)) sources <- c(sources, transcriptomics = "Transcriptomics")
      if (!is.null(private$dde_prot)) sources <- c(sources, proteomics = "Proteomics")
      if (!is.null(private$dde_metabo)) sources <- c(sources, metabolomics = "Metabolomics")
      return(sources)
    }
  )
)

MovidaModelSingleton <- local({
  instance <- NULL

  list(
    get_instance = function(movida_list = NULL) {
      if (is.null(instance)) {
        if (is.null(movida_list)) {
          stop("First call must provide movida_list to initialize the singleton.")
        }
        instance <<- MovidaModel$new(movida_list)
      }
      instance
    }
  )
})
