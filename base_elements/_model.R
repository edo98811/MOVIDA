library(R6)

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
    metadata = NULL
  ),
  
  public = list(
    initialize = function(movida_list) {
      private$annotation_df_prot <- movida_list$annotation_df_prot
      private$annotation_df_trans <- movida_list$annotation_df_trans
      private$annotation_df_metabo <- movida_list$annotation_df_metabo
      private$results_prot <- movida_list$results_prot
      private$results_trans <- movida_list$results_trans
      private$results_metabo <- movida_list$results_metabo
      private$se_prot <- movida_list$se_prot
      private$se_trans <- movida_list$se_trans
      private$se_metabo <- movida_list$se_metabo
      private$metadata <- movida_list$metadata
    },
    
    get_de = function(type = c("prot", "trans", "metabo")) {
      type <- match.arg(type)
      # Example: return differential expression results based on type
      if (type == "prot") return(private$results_prot)
      if (type == "trans") return(private$results_trans)
      if (type == "metabo") return(private$results_metabo)
    },
    
    get_enrichment = function() {
      # Placeholder: implement enrichment logic here
      message("Running enrichment analysis...")
      # Return some result or data frame
    },
    
    get_expression = function(type = c("prot", "trans", "metabo")) {
      type <- match.arg(type)
      # Placeholder: return expression data, e.g. the se_... slots
      if (type == "prot") return(private$se_prot)
      if (type == "trans") return(private$se_trans)
      if (type == "metabo") return(private$se_metabo)
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