#' @title MovidaModel Class
#' @description An R6 class for managing and analyzing multi-omics data.
#' @details The `MovidaModel` class provides methods for initializing, processing,
#' and retrieving relationships and data from proteomics, transcriptomics,
#' and metabolomics experiments. It supports parallel processing and
#' includes functionality for filtering and retrieving data based on
#' specific criteria.
#'
#' @section Private Fields:
#' \describe{
#'   \item{results_prot}{Results from proteomics analysis.}
#'   \item{results_trans}{Results from transcriptomics analysis.}
#'   \item{results_metabo}{Results from metabolomics analysis.}
#'   \item{se_prot}{SummarizedExperiment object for proteomics data.}
#'   \item{se_trans}{SummarizedExperiment object for transcriptomics data.}
#'   \item{se_metabo}{SummarizedExperiment object for metabolomics data.}
#'   \item{contrasts}{List of contrasts used in the analysis.}
#'   \item{inchi_to_ensembl}{Relationships between InChI and Ensembl IDs.}
#'   \item{inchi_to_uniprot}{Relationships between InChI and UniProt IDs.}
#'   \item{uniprot_to_ensembl}{Relationships between UniProt and Ensembl IDs.}
#'   \item{organism}{Organism information for the data.}
#'   \item{determine_contrasts}{Function to determine contrasts from the input data.}
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
#' - `get_dea(target, contrast, FDRpvalue = NULL, FDRadj = NULL)`: Retrieves differential expression results
#'   for the given data target and contrast, with optional FDR-based filtering.
#'
#' - `get_fea(target, contrast, FDRpvalue = NULL, FDRadj = NULL)`: Retrieves functional enrichment analysis results
#'   for the specified target and contrast, with optional filtering.
#'
#' - `get_values_all(target, return_se = FALSE)`: Returns all assay data (or `SummarizedExperiment` object)
#'   for the selected data target.
#'
#' - `get_values_features(features, source = NULL)`: Retrieves expression data for specific features.
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
#'   se_prot = NULL,
#'   se_trans = NULL,
#'   se_metabo = NULL,
#'   organism = "Hs"
#' )
#' model <- MovidaModel$new(movida_list)
#'
library(R6)
library(SummarizedExperiment)
source("utils/utils_find_relationships.R")
source("utils/utils_validate_input.R")

MovidaModel <- R6Class("MovidaModel",
  private = list(
    results_prot = NULL,
    results_trans = NULL,
    results_metabo = NULL,
    se_prot = NULL,
    se_trans = NULL,
    se_metabo = NULL,
    contrasts = NULL,
    inchi_to_ensembl = NULL,
    inchi_to_uniprot = NULL,
    uniprot_to_ensembl = NULL,
    organism = NULL,
    determine_contrasts = function(movida_list) {
      if (!"contrasts" %in% names(movida_list) || is.null(movida_list$contrasts)) {
        private$contrasts <- union(
          names(movida_list$results_prot),
          union(
            names(movida_list$results_trans),
            names(movida_list$results_metabo)
          )
        )
        warning("Contrasts not provided in movida_list. Using union of names from results_prot, results_trans, and results_metabo.")
      } else {
        private$contrasts <- movida_list$contrasts
      }
    }
  ),
  public = list(
    initialize = function(movida_list) {
      check_movida_list(movida_list)
      # Initialize private variables with data from the provided movida_list
      private$results_prot <- movida_list$results_prot
      private$results_trans <- movida_list$results_trans
      private$results_metabo <- movida_list$results_metabo
      private$se_prot <- movida_list$se_prot
      private$se_trans <- movida_list$se_trans
      private$se_metabo <- movida_list$se_metabo
      private$organism <- movida_list$organism
      # self$build_relationships()
      private$determine_contrasts(movida_list)
    },
    build_relationships = function() {
      # Use parallel processing for building relationships!! used for inchi relationships
      if (inherits(future::plan(), "multisession")) {
        future::plan(future::sequential) # Stop the existing multisession plan
        warning("Parallel processing was already set to multisession. Restarting it.")
      }
      future::plan(future::multisession, workers = parallel::detectCores(logical = FALSE)) # Set the number of workers as needed

      # Build inchi relationships if metabolomics data is available
      if (!is.null(private$se_metabo)) {
        # Build inchi to Ensembl relationships
        message("Creating relationships dataset inchi_to_uniprot...")
        if (!is.null(private$se_prot)) {
          private$inchi_to_ensembl <- build_inchi_relationships(rowData(private$se_metabo), rowData(private$se_prot), get_ensembl = FALSE, get_uniprot = TRUE)
        }
        message("Done!")

        # Build inchi to UniProt relationships
        message("Creating relationships dataset inchi_to_ensembl...")
        if (!is.null(private$se_trans)) {
          private$inchi_to_uniprot <- build_inchi_relationships(rowData(private$se_metabo), rowData(private$se_trans), get_ensembl = TRUE, get_uniprot = FALSE)
        }
        message("Done!")
      }

      # Build UniProt-to-Ensembl relationships if both proteomics and transcriptomics data are available
      if (!is.null(private$se_prot) && !is.null(private$se_trans)) {
        private$uniprot_to_ensembl <- build_uniprot_to_ensembl(rowData(private$se_prot), rowData(private$se_trans), private$organism)
      }
    },
    load_relationships = function(folder_path) {
      if (!is.character(folder_path) || length(folder_path) != 1) {
        stop("Argument 'folder_path' must be a single string.")
      }
      if (!dir.exists(folder_path)) {
        stop("The specified folder does not exist.")
      }

      private$inchi_to_ensembl <- readRDS(file.path(folder_path, "inchi_to_ensembl.rds"))
      private$inchi_to_uniprot <- readRDS(file.path(folder_path, "inchi_to_uniprot.rds"))
      private$uniprot_to_ensembl <- readRDS(file.path(folder_path, "uniprot_to_ensembl.rds"))
      return(invisible(TRUE))
    },
    save_relationships = function(folder_path) {
      if (!is.character(folder_path) || length(folder_path) != 1) {
        stop("Argument 'folder_path' must be a single string.")
      }
      if (!dir.exists(folder_path)) {
        dir.create(folder_path, recursive = TRUE)
      }

      saveRDS(private$inchi_to_ensembl, file.path(folder_path, "inchi_to_ensembl.rds"))
      saveRDS(private$inchi_to_uniprot, file.path(folder_path, "inchi_to_uniprot.rds"))
      saveRDS(private$uniprot_to_ensembl, file.path(folder_path, "uniprot_to_ensembl.rds"))
    },
    get_contrasts = function() {
      return(private$contrasts)
    },
    get_inchi_to_ensembl = function() {
      return(private$inchi_to_ensembl)
    },
    get_inchi_to_uniprot = function() {
      return(private$inchi_to_uniprot)
    },
    get_uniprot_to_ensembl = function() {
      return(private$uniprot_to_ensembl)
    },
    get_dea = function(source, contrast, FDRpvalue = NULL, FDRadj = NULL) {
      if (!is.character(source) || length(source) != 1) {
        stop("Argument 'source' must be a single string.")
      }
      if (!is.character(contrast) || length(contrast) != 1) {
        stop("Argument 'contrast' must be a single string.")
      }

      if (source == "proteomics") {
        data <- private$results_prot[[contrast]]$tbl_res_all
      } else if (source == "transcriptomics") {
        data <- private$results_trans[[contrast]]$tbl_res_all
      } else if (source == "metabolomics") {
        data <- private$results_metabo[[contrast]]$tbl_res_all
      } else {
        stop("Invalid source. Must be one of 'proteomics', 'transcriptomics', or 'metabolomics'.")
      }

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
    get_fea = function(source, contrast, FDRpvalue = NULL, FDRadj = NULL) {
      if (!is.character(source) || length(source) != 1) {
        stop("Argument 'source' must be a single string.")
      }
      if (!is.character(contrast) || length(contrast) != 1) {
        stop("Argument 'contrast' must be a single string.")
      }

      if (source == "proteomics") {
        data <- private$results_prot[[contrast]]$topGO_tbl
      } else if (source == "transcriptomics") {
        data <- private$results_trans[[contrast]]$topGO_tbl
      } else if (source == "metabolomics") {
        data <- private$results_metabo[[contrast]]$topGO_tbl
      } else {
        stop("Invalid source. Must be one of 'proteomics', 'transcriptomics', or 'metabolomics'.")
      }

      if (!is.null(FDRpvalue)) {
        data <- data[data$FDRpvalue <= FDRpvalue, ]
      }
      if (!is.null(FDRadj)) {
        data <- data[data$FDRadj <= FDRadj, ]
        if (!is.null(FDRpvalue)) {
          warning("Both FDRpvalue and FDRadj are set. Filtering will be applied based on FDRadj")
        }
      }
      if (!is.data.frame(data)) {
        stop("The filtered data is not a data frame. Please check the input data structure.")
      }
      return(data)
    },
    get_values_all = function(source, return_se = FALSE) {
      # Validate that 'source' is a single string
      if (!is.character(source) || length(source) != 1) {
        stop("Argument 'source' must be a single string.")
      }

      # Select the appropriate SummarizedExperiment object based on the source
      if (source == "proteomics") {
        se <- private$se_prot
      } else if (source == "transcriptomics") {
        se <- private$se_trans
      } else if (source == "metabolomics") {
        se <- private$se_metabo
      } else {
        stop("Invalid source. Must be one of 'proteomics', 'transcriptomics', or 'metabolomics'.")
      }

      # Return the assay data or the SummarizedExperiment object based on 'return_se'
      if (return_se) {
        return(se) # Return the assay data (matrix of values)
      } else {
        return(assay(se)) # Return the SummarizedExperiment object
      }
    },
    get_related_features = function(feature, target) {
      if (!is.character(c(feature))) {
        stop("Argument 'features' must be a vector of strings.")
      }

      # Validate that the features are of a recognized type (Ensembl, inchi, or UniProt)
      if (!suppressWarnings(check_ensembl(c(feature))) && !suppressWarnings(check_inchi(c(feature))) && !suppressWarnings(check_uniprot(c(feature)))) {
        stop("Error: rownames(rowData(se)) must be of type Ensembl, inchi, or UniProt.")
      }

      # Check if 'target' is a single string
      if (!is.character(target) || length(target) != 1) {
        stop("Argument 'target' must be a single string.")
      }

      # Return the related features based on the target type
      if (suppressWarnings(check_ensembl(c(feature)))) {
        if (target == "proteomics") {
          return(private$uniprot_to_ensembl[private$uniprot_to_ensembl$ENSEMBL %in% feature, , drop = FALSE]$UNIPROT)
        } else if (target == "transcriptomics") {
          return(c(feature)) # Ensembl to Ensembl is a direct match
        } else if (target == "metabolomics") {
          return(private$inchi_to_ensembl[private$inchi_to_ensembl$ENSEMBL %in% feature, , drop = FALSE]$INCHIKEY)
        }
      } else if (suppressWarnings(check_inchi(c(feature)))) {
        if (target == "proteomics") {
          return(private$inchi_to_uniprot[private$inchi_to_uniprot$INCHIKEY %in% feature, , drop = FALSE]$UNIPROT)
        } else if (target == "transcriptomics") {
          return(private$inchi_to_ensembl[private$inchi_to_ensembl$INCHIKEY %in% feature, , drop = FALSE]$ENSEMBL)
        } else if (target == "metabolomics") {
          return(c(feature)) # inchi to inchi is a direct match
        }
      } else if (suppressWarnings(check_uniprot(c(feature)))) {
        if (target == "proteomics") {
          return(c(feature)) # UniProt to UniProt is a direct match
        } else if (target == "transcriptomics") {
          return(private$uniprot_to_ensembl[private$uniprot_to_ensembl$UNIPROT %in% feature, , drop = FALSE]$ENSEMBL)
        } else if (target == "metabolomics") {
          return(private$inchi_to_uniprot[private$inchi_to_uniprot$UNIPROT %in% feature, , drop = FALSE]$INCHIKEY)
        }
      }
    },
    get_values_features = function(features, source = NULL) {
      # Check if 'features' is a vector of strings
      if (!is.character(features) || !is.vector(features)) {
        stop("Argument 'features' must be a vector of strings.")
      }

      # Determine the source type based on the input features
      if (check_ensembl(features)) {
        source <- "transcriptomics"
      } else if (check_inchi(features)) {
        source <- "metabolomics"
      } else if (check_uniprot(features)) {
        source <- "proteomics"
      } else {
        stop("Error: features must be of type Ensembl, inchi, or UniProt.")
      }

      # Select the appropriate SummarizedExperiment object based on the source
      if (source == "proteomics") {
        se <- private$se_prot
      } else if (source == "transcriptomics") {
        se <- private$se_trans
      } else if (source == "metabolomics") {
        se <- private$se_metabo
      } else {
        stop("Invalid source type. Must be one of 'proteomics', 'transcriptomics', or 'metabolomics'.")
      }

      # Return the assay data for the specified features
      return(assay(se)[features, , drop = FALSE])
    },
    get_values_group = function(groups, source) {
      # Check if 'groups' is a vector of strings
      if (!is.character(groups) || !is.vector(groups)) {
        stop("Argument 'groups' must be a vector of strings.")
      }

      # Select the appropriate SummarizedExperiment object based on the source
      if (source == "proteomics") {
        se_source <- private$se_prot
      } else if (source == "transcriptomics") {
        se_source <- private$se_trans
      } else if (source == "metabolomics") {
        se_source <- private$se_metabo
      } else {
        stop("Invalid source type. Must be one of 'proteomics', 'transcriptomics', or 'metabolomics'.")
      }

      # Check that all groups in 'groups' are present in colData(se_source)$group
      if (!all(groups %in% colData(se_source)$group)) {
        stop("Some groups in 'groups' are not present in colData(se_source)$group.")
      }

      # Subset se_source based on the specified groups
      subset_indices <- colData(se_source)$group %in% groups
      return(se_source[, subset_indices])
    },
    get_values_samples = function(samples, source) {
      # Check if 'samples' is a vector of strings
      if (!is.character(samples) || !is.vector(samples)) {
        stop("Argument 'samples' must be a vector of strings.")
      }

      # Select the appropriate SummarizedExperiment object based on the source
      if (source == "proteomics") {
        se_source <- private$se_prot
      } else if (source == "transcriptomics") {
        se_source <- private$se_trans
      } else if (source == "metabolomics") {
        se_source <- private$se_metabo
      } else {
        stop("Invalid source type. Must be one of 'proteomics', 'transcriptomics', or 'metabolomics'.")
      }

      # Check that all samples in 'samples' are present in colnames(se_source)
      if (!all(samples %in% colnames(se_source))) {
        stop("Some samples in 'samples' are not present in colnames(se_source).")
      }

      # Subset se_source based on the specified samples
      subset_indices <- colnames(se_source) %in% samples
      return(se_source[, subset_indices])
    },
    get_features_list = function(source) {
      if (source == "proteomics") {
        return(as.data.frame(rownames(private$se_prot)))
      } else if (source == "transcriptomics") {
        return(as.data.frame(rownames(private$se_trans)))
      } else if (source == "metabolomics") {
        return(as.data.frame(rownames(private$se_metabo)))
      } else {
        stop("Invalid source. Must be one of 'proteomics', 'transcriptomics', or 'metabolomics'.")
      }
    },
    get_metadata_columns = function(source) {
      if (source == "proteomics") {
        return(colnames(private$se_prot))
      } else if (source == "transcriptomics") {
        return(colnames(private$se_trans))
      } else if (source == "metabolomics") {
        return(colnames(private$se_metabo))
      } else {
        stop("Invalid source. Must be one of 'proteomics', 'transcriptomics', or 'metabolomics'.")
      }
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
