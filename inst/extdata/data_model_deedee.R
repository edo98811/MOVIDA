
MovidaModel <- R6Class("MovidaModel",
  private = list(
    annotation_df = NULL,
    annotation_df_metabo = NULL,
    DeeDee_prot = NULL,
    DeeDee_trans = NULL,
    DeeDee_metabo = NULL,
    metadata = NULL,
    contrasts = NULL
  ),

  public = list(
    initialize = function(movida_list) {
      if (!"contrasts" %in% names(movida_list)) {
      private$contrasts <- intersect(          
          names(movida_list$results_prot), 
          names(movida_list$results_trans)
        # intersect(Reduce(intersect, list(
        #   names(movida_list$results_prot), 
        #   names(movida_list$results_trans), 
        #   names(movida_list$results_metabo)
        # ))
      )
        warning("Contrasts not provided in movida_list. Using intersection of names from results_prot, results_trans, and results_metabo.")
      } else {
        private$contrasts <- movida_list$contrasts
      }
    if (is.null(movida_list$annotation_df)) {
      private$annotation_df <- create_annotation(movida_list)
      warning("annotation_df not provided in movida_list. Created using create_annotation function.")
    } else {
      private$annotation_df <- movida_list$annotation_df
    }
    private$annotation_df_metabo <- movida_list$annotation_df_metabo,
    private$DeeDee_prot <- movida_list$DeeDee_prot,
    private$DeeDee_trans <- movida_list$DeeDee_trans,
    private$DeeDee_metabo <- movida_list$DeeDee_metabo,
    private$metadata <- movida_list$metadata,
    },
    get_contrasts = function() {
      return(private$contrasts)
    },
    get_de = function(type, contrast, FDRpvalue = NULL, FDRadj = NULL) {
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

    get_enrichment = function(type, contrast, FDRpvalue = NULL, FDRadj = NULL) {
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

    get_expression = function(type, return_se = FALSE) {

      if (!is.character(type) || length(type) != 1) {
        stop("Argument 'type' must be a single string.")
      }

      if (type == "proteomics") {
        se <- private$se_prot
      } else if (type == "transcriptomics") {
        se <- private$se_trans
      } else if (type == "metabolomics") {
        se <- private$se_metabo
      } else {
        stop("Invalid type. Must be one of 'prot', 'trans', or 'metabo'.")
      }

      if (return_se) return(assay(se))
      else return(se)
    },
    get_anno_df = function(type) {

      if (!is.character(type) || length(type) != 1) {
        stop("Argument 'type' must be a single string.")
      }

      if (type == "proteomics") {
        anno_df <- private$annotation_df_prot
      } else if (type == "transcriptomics") {
        anno_df <- private$annotation_df_trans
      } else if (type == "metabolomics") {
        anno_df <- private$annotation_df_metabo
      } else {
        stop("Invalid type. Must be one of 'prot', 'trans', or 'metabo'.")
      }

      return(anno_df)
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