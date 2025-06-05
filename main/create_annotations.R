#' Create Annotations for SummarizedExperiment Object
#'
#' This function generates annotations for a given `SummarizedExperiment` object 
#' by mapping identifiers between different annotation types (e.g., SYMBOL, ENSEMBL, ENTREZID, UNIPROT).
#'
#' @param params A list containing parameters for the function. Must include the `species` key, 
#'        which specifies the species. Supported values are `"Mm"` for mouse and `"Hs"` for human.
#' @param se A `SummarizedExperiment` object. The row names of this object are used as keys for annotation mapping.
#' @param source_type A character string specifying the type of identifier used in the row names of `se`. 
#'        Supported values include `"ENSEMBL"`, `"SYMBOL"`, `"UNIPROT"`, etc.
#' @param columns A character vector specifying the types of annotations to retrieve. 
#'        Default is `c("SYMBOL", "ENSEMBL", "ENTREZID", "UNIPROT")`.
#' @param force_creation A logical value indicating whether to force the creation of annotations 
#'        even if they already exist in the parent environment. Default is `FALSE`.
#'
#' @return A data frame containing the requested annotations. The row names of the returned data frame 
#'         match the row names of the input `SummarizedExperiment` object, and the column names correspond 
#'         to the requested annotation types.
#'
#' @details
#' The function uses the `AnnotationDbi` package to map identifiers between different annotation types. 
#' If any entries cannot be mapped, a warning is issued, and `NA` values are returned for those entries. 
#' If all entries fail to map, the function stops with an error.
#'
#' @examples
#' \dontrun{
#' library(SummarizedExperiment)
#' library(org.Mm.eg.db)
#' 
#' # Example SummarizedExperiment object
#' se <- SummarizedExperiment(assays = list(counts = matrix(1:4, nrow = 2)),
#'                             rowData = DataFrame(ENSEMBL = c("ENSMUSG00000000001", "ENSMUSG00000000003")))
#' 
#' # Create annotations
#' params <- list(species = "Mm")
#' annotations <- create_annotations(params, se, source_type = "ENSEMBL", 
#'                                    columns = c("SYMBOL", "ENTREZID"))
#' print(annotations)
#' }
#'
#' @importFrom AnnotationDbi mapIds keytypes
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @export

create_annotations <- function(params, se, source_type = "ENSEMBL", columns = c("SYMBOL", "ENSEMBL", "ENTREZID", "UNIPROT"), force_creation = FALSE) {

  # Validate the input parameters
  if (!is.list(params) || !all(c("species") %in% names(params))) {
    stop("Invalid parameters provided. Please provide a list with 'species', 'analysis_folder'.")
  }
  species <- params$species

  # Check if annotations already exist in the parent environment and return them if force_creation is FALSE
  if (exists("anns", envir = parent.frame()) && !force_creation) {
    message("Annotations already exist in the parent environment. Returning existing annotations...")
    return(get("anns", envir = parent.frame()))
  } 

  # Validate that the input object is a SummarizedExperiment object
  if (!inherits(se, "SummarizedExperiment")) {
    stop("The provided object is not a SummarizedExperiment object.")
  }

  # Select the appropriate organism database based on the species parameter
  if (species == "Mm") {
    orgdb <- org.Mm.eg.db::org.Mm.eg.db  # Use mouse-specific annotation database
  } else if (species == "Hs") {
    orgdb <- org.Hs.eg.db::org.Hs.eg.db  # Use human-specific annotation database
  } else {
    stop("Unsupported species. Please use 'Mm' for mouse or 'Hs' for human.")  # Handle unsupported species
  }

  # Validate source_type and columns against available keys/columns in orgdb
  valid_keys <- AnnotationDbi::keytypes(orgdb)

  if (!(source_type %in% valid_keys)) {
    stop(sprintf("source_type '%s' is not a valid keytype for the selected species/orgdb. Valid options include: %s", 
                 source_type, paste(valid_keys, collapse = ", ")))
  }
  
  # Check requested columns against keytypes as well
  invalid_cols <- setdiff(columns, valid_keys)
  if (length(invalid_cols) > 0) {
    stop(sprintf("The following requested columns are invalid for the selected orgdb: %s", paste(invalid_cols, collapse = ", ")))
  }
  #  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
  #  [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
  # [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MGI"         
  # [16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
  # [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UNIPROT"   

  # Map each requested annotation column
  anns <- do.call(
    cbind.data.frame,
    lapply(columns, function(dest_col) {
      create_column(keys_list, source_type, dest_col, orgdb)
    }),
    stringsAsFactors = FALSE
  )

  # to preserve the names of the rows in the original SummarizedExperiment object and have tehe right colnames
  colnames(anns) <- columns
  rownames(anns) <- rownames(se)

  return(anns)
}

create_column <- function(keys_list, source_type, destination_type, orgdb) {

  # Create a new column in the data frame based on the keys_list
  if(source_type == destination_type) {
    return(keys_list)
  }

  # Create a new column in the data frame based on the keys_list
  new_column <- AnnotationDbi::select(
        orgdb,
        keys = keys_list,
        columns = destination_type,
        keytype = source_type
      )

  # Check for missing values and warn if any entries were not found
  num_missing <- sum(is.na(new_column))
  if (num_missing > 0) {
    warning(sprintf("Mapping from %s to %s: %d/%d entries not found (returned NA).",
                    source_type, destination_type, num_missing, length(keys_list)))
  }

  # If all entries are missing, stop the function with an error
  if (num_missing == length(keys_list)) {
    stop(sprintf("Mapping from %s to %s failed: No entries could be mapped. Did you set the correct parameters? ", source_type, destination_type))
  }

  return(new_column)
}