library(R6)
source("utils/create_annotations.R")
source("utils/validate_input.R")

MovidaModel <- R6Class("MovidaModel",
  private = list(
    annotation_df_prot = NULL,
    annotation_df_trans = NULL,
    annotation_df_metabo = NULL,
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
      if (!"contrasts" %in% names(movida_list)) {
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
      private$organism <- if (check_organism(movida_list$organism)) movida_list$organism else stop("Invalid organism provided in movida_list.")
      self$build_relationships()
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
        if (!is.null(private$se_prot)) {
          private$inchi_to_ensembl <- inchi_relationships(
            private$se_metabo, private$se_prot,
            get_ensembl = FALSE, get_uniprot = TRUE
          )
        }
        # Build inchi to UniProt relationships
        if (!is.null(private$se_trans)) {
          private$inchi_to_uniprot <- inchi_relationships(
            private$se_metabo, private$se_prot,
            get_ensembl = TRUE, get_uniprot = FALSE
          )
        }
      }

      # Build UniProt-to-Ensembl relationships if both proteomics and transcriptomics data are available
      if (!is.null(private$se_prot) && !is.null(private$se_trans)) {
        private$uniprot_to_ensembl <- uniprot_relationships(
          private$se_prot, private$se_trans,
          organism = private$organism
        )
      }
    },
    load_relationships = function(folder_path) {
      if (!dir.exists(folder_path)) {
      stop("The specified folder does not exist.")
      }
      
      private$inchi_to_ensembl <- readRDS(file.path(folder_path, "inchi_to_ensembl.rds"))
      private$inchi_to_uniprot <- readRDS(file.path(folder_path, "inchi_to_uniprot.rds"))
      private$uniprot_to_ensembl <- readRDS(file.path(folder_path, "uniprot_to_ensembl.rds"))
    },
    save_relationships = function(folder_path) {
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
    get_dea = function(type, contrast, FDRpvalue = NULL, FDRadj = NULL) {
      if (!is.character(type) || length(type) != 1) {
        stop("Argument 'type' must be a single string.")
      }
      if (!is.character(contrast) || length(contrast) != 1) {
        stop("Argument 'contrast' must be a single string.")
      }

      if (type == "proteomics") {
        data <- private$results_prot[[contrast]]$tbl_res_all
      } else if (type == "transcriptomics") {
        data <- private$results_trans[[contrast]]$tbl_res_all
      } else if (type == "metabolomics") {
        data <- private$results_metabo[[contrast]]$tbl_res_all
      } else {
        stop("Invalid type. Must be one of 'prot', 'trans', or 'metabo'.")
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
    get_fea = function(type, contrast, FDRpvalue = NULL, FDRadj = NULL) {
      if (!is.character(type) || length(type) != 1) {
        stop("Argument 'type' must be a single string.")
      }
      if (!is.character(contrast) || length(contrast) != 1) {
        stop("Argument 'contrast' must be a single string.")
      }

      if (type == "proteomics") {
        data <- private$results_prot[[contrast]]$topGO_tbl
      } else if (type == "transcriptomics") {
        data <- private$results_trans[[contrast]]$topGO_tbl
      } else if (type == "metabolomics") {
        data <- private$results_metabo[[contrast]]$topGO_tbl
      } else {
        stop("Invalid type. Must be one of 'prot', 'trans', or 'metabo'.")
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
    get_values_all = function(type, return_se = FALSE) {
      # Validate that 'type' is a single string
      if (!is.character(type) || length(type) != 1) {
        stop("Argument 'type' must be a single string.")
      }

      # Select the appropriate SummarizedExperiment object based on the type
      if (type == "proteomics") {
        se <- private$se_prot
      } else if (type == "transcriptomics") {
        se <- private$se_trans
      } else if (type == "metabolomics") {
        se <- private$se_metabo
      } else {
        stop("Invalid type. Must be one of 'prot', 'trans', or 'metabo'.")
      }

      # Return the assay data or the SummarizedExperiment object based on 'return_se'
      if (return_se) {
        return(assay(se)) # Return the assay data (matrix of values)
      } else {
        return(se) # Return the SummarizedExperiment object
      }
    },
    get_related_features = function(feature, target) {
      if (!is.character(c(feature))) {
        stop("Argument 'features' must be a vector of strings.")
      }

      # Validate that the features are of a recognized type (Ensembl, inchi, or UniProt)
      if (!check_ensembl(c(feature)) && !check_inchi(c(feature)) && !check_uniprot(c(feature))) {
        stop("Error: rownames(rowData(se)) must be of type Ensembl, inchi, or UniProt.")
      }

      # Check if 'type' is a single string
      if (!is.character(type) || length(type) != 1) {
        stop("Argument 'type' must be a single string.")
      }

      # Return the related features based on the target type
      if (check_ensembl(c(feature))) {
        if (target == "uniprot") {
          return(private$uniprot_to_ensembl_dataframe[private$uniprot_to_ensembl_dataframe$ENSEMBL %in% feature, , drop = FALSE]$UNIPROT)
        } else if (target == "ensembl") {
          return(feature) # Ensembl to Ensembl is a direct match
        } else if (target == "inchi") {
          return(private$inchi_to_ensembl[private$inchi_to_ensembl$ENSEMBL %in% feature, , drop = FALSE]$inchi)
        }
      } else if (check_inchi(c(feature))) {
        if (target == "uniprot") {
          return(private$inchi_to_uniprot[private$inchi_to_uniprot$inchi %in% feature, , drop = FALSE]$UNIPROT)
        } else if (target == "ensembl") {
          return(private$inchi_to_ensembl[private$inchi_to_ensembl$inchi %in% feature, , drop = FALSE]$ENSEMBL)
        } else if (target == "inchi") {
          return(feature) # inchi to inchi is a direct match
        }
      } else if (check_uniprot(c(feature))) {
        if (target == "uniprot") {
          return(feature) # UniProt to UniProt is a direct match
        } else if (target == "ensembl") {
          return(private$uniprot_to_ensembl_dataframe[private$uniprot_to_ensembl_dataframe$UNIPROT %in% feature, , drop = FALSE]$ENSEMBL)
        } else if (target == "inchi") {
          return(private$inchi_to_uniprot[private$inchi_to_uniprot$UNIPROT %in% feature, , drop = FALSE]$inchi)
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
