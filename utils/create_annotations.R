# Import relevant libraries
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(httr2)
library(jsonlite)
library(future)
library(furrr)

# List of functions in this file:
# 1. uniprot_chebi_query_async: Asynchronously queries UniProt for multiple ChEBI IDs and filters results based on target data.
# 2. uniprot_chebi_query: Queries UniProt for a single ChEBI ID and retrieves associated protein and gene IDs.
# 3. build_chebi_relationships: Builds relationships between ChEBI IDs, UniProt IDs, and Ensembl IDs based on metabolomics and target data.
# 4. build_uniprot_to_ensembl: Maps UniProt IDs to Ensembl IDs for a given organism and filters results based on transcript data.
# 5. add_chebi_ensembl_relationships: NOT USED Merges ChEBI-to-UniProt relationships with UniProt-to-Ensembl mappings to establish ChEBI-to-Ensembl relationships. 
# 6. chebi_relationships: Generates ChEBI relationships using SummarizedExperiment objects for metabolomics and target data.
# 7. uniprot_relationships: Builds UniProt-to-Ensembl relationships using SummarizedExperiment objects for protein and transcript data.

# Wrapper for async querying multiple ChEBI IDs
uniprot_chebi_query_async <- function(chebi_ids, rowdata_target, get_ensembl = FALSE, get_uniprot = TRUE, sleep_time = 0.05) {
  results <- furrr::future_map(chebi_ids, function(id) { # like purrr map but works in parallel: https://www.rdocumentation.org/packages/furrr/versions/0.3.1/topics/future_map
    Sys.sleep(sleep_time) # Avoid hitting the API too hard (needed?)
    tryCatch(
      {
        # Query UniProt for each ChEBI ID (fucntion return list)
        query_result <- uniprot_chebi_query(id, get_ensembl = get_ensembl, get_uniprot = get_uniprot)

        # Filter the UniProt IDs to include only those present in the rownames of rowdata_target
        if (get_uniprot && !is.null(query_result$protein_ids)) {
          query_result$protein_ids <- query_result$protein_ids[query_result$protein_ids %in% rownames(rowdata_target)]
        }
        # Filter the Ensembl IDs to include only those present in the rownames of rowdata_target
        if (get_ensembl && !is.null(query_result$gene_ids)) {
          query_result$gene_ids <- query_result$gene_ids[query_result$gene_ids %in% rownames(rowdata_target)]
        }

        query_result
      },
      error = function(e) {
        warning(sprintf("Failed for %s: %s", id, e$message))
        return(list(protein_ids = NULL, gene_ids = NULL))
      }
    )
  }, .options = furrr_options(seed = TRUE)) # Ensures reproducibility

  names(results) <- chebi_ids
  return(results)
}

#' @return A list containing the following elements:
#' \describe{
#'   \item{protein_ids}{A character vector of UniProt protein IDs extracted from the API response.
#'   If \code{get_uniprot} is \code{FALSE}, this will be \code{NULL}.}
#'   \item{gene_ids}{A character vector of Ensembl gene IDs extracted from the API response.
#'   If \code{get_ensembl} is \code{FALSE}, this will be \code{NULL}.
#'   If no Ensembl references are found for a result, the corresponding value will be \code{NA}.}
#' }
# Function to query UniProt for ChEBI-related data
uniprot_chebi_query <- function(chebi_id, get_ensembl = FALSE, get_uniprot = TRUE) {
  if (!check_chebi(chebi_id)) {
    stop("Not a valid ChEBI ID: ", chebi_id)
  }
  # Base URL for the UniProt REST API
  base_url <- "https://rest.uniprot.org/uniprotkb/search"

  # Define the query parameters for the API request


  params <- list(
    query = paste0(chebi_id, " AND reviewed:true"), # Example query; replace with actual ChEBI ID logic
    fields = paste(c(
      if (get_uniprot) "primaryAccession",
      if (get_ensembl) "xref_ensembl"
    ), collapse = ","), # Fields to retrieve
    sort = "accession desc" # Sort results by accession in descending order
  )

  # Create the request object
  req <- request(base_url)
  req <- req |> req_headers(accept = "application/json")
  req <- req |> req_url_query(!!!params)

  resp <- req_perform(req) # Perform the API request

  # Handle non-200 response
  status_code <- resp_status(resp)
  if (status_code != 200) {
    error_body <- tryCatch(
      jsonlite::fromJSON(resp_body_string(resp))$messages,
      error = function(e) substr(resp_body_string(resp), 1, 300) # fallback
    )
    stop(sprintf("UniProt API error (HTTP %d): %s", status_code, paste(error_body, collapse = "; ")))
  }

  # Parse the JSON response into an R object
  data <- fromJSON(content(resp, "text", encoding = "UTF-8"))
  if (length(data$results) == 0) {
    warning("No matches found for the provided ChEBI ID.")
    return(list(protein_ids = NULL, gene_ids = NULL))
    # stop("No matches found for the provided ChEBI ID during mapping to UniProt.")
  }

  # Extract UniProt protein IDs if requested  # Extract Ensembl gene IDs if requested
  # STRUCUTRE OF THE UNIPROT REPLY:
  #   {
  #   "entryType": "UniProtKB reviewed (Swiss-Prot)",
  #   "primaryAccession": "P50207",
  #   "uniProtKBCrossReferences": [
  #     {
  #       "database": "Ensembl",
  #       "id": "ENSMUST00000001700.7",
  #       "properties": [
  #         {
  #           "key": "ProteinId",
  #           "value": "ENSMUSP00000001700.7"
  #         },
  #         {
  #           "key": "GeneId",
  #           "value": "ENSMUSG00000001655.7"
  #         }
  #       ]
  #     }
  #   ],
  #   "extraAttributes": {
  #     "uniParcId": "UPI0000020B9F"
  #   }
  # },

  if (get_uniprot) {
    protein_ids <- sapply(data$results, function(x) {
      if (is.null(x$primaryAccession)) { # If primaryAccession is NULL, return NA
        return(NA)
      } else {
        x$primaryAccession
      }
    })
  } else {
    protein_ids <- NULL
  }

  # Using Filter: https://rdrr.io/r/base/funprog.html. extracts the elements of a vector for which a predicate (logical) function gives true.
  if (get_ensembl) {
    gene_ids <- lapply(data$results, function(x) {
      if (is.null(x$uniProtKBCrossReferences)) {
        return(NA)
      } # If uniProtKBCrossReferences is NULL, return NA
      ensembl_refs <- Filter(function(ref) ref$database == "Ensembl", x$uniProtKBCrossReferences)
      gene_ids <- unlist(lapply(ensembl_refs, function(ref) {
        props <- ref$properties # Extract properties from the Ensembl cross-reference as list
        genes <- Filter(function(p) p$key == "GeneId", props) # If the key is "GeneId", then save it in genes
        sapply(genes, function(g) sub("\\.\\d+$", "", g$value)) # Return the value of all the GeneId properties, removing version numbers
      }))
      if (length(gene_ids) == 0) { # If no Ensembl references are found, return NA
        return(NA)
      }
      return(unique(gene_ids))
    })
  } else {
    gene_ids <- NULL
  }

  # Return a list containing the extracted protein and gene IDs
  return(list(protein_ids = protein_ids, gene_ids = gene_ids))
}

#' Build ChEBI Relationships
#'
#' This function generates a data frame of relationships between ChEBI IDs, UniProt IDs,
#' and Ensembl IDs based on metabolomics and target data. It queries UniProt to retrieve
#' protein and gene information associated with the provided ChEBI IDs and filters the
#' results based on their presence in the target data.
#'
#' @param rowdata_metabo A data frame or similar object containing metabolomics data.
#'        This should include ChEBI IDs.
#' @param rowdata_target A data frame or similar object containing target data.
#'        The row names should correspond to UniProt or Ensembl IDs.
#' @param get_ensembl Logical. If TRUE, Ensembl gene IDs will be retrieved and included
#'        in the output. Default is TRUE.
#' @param get_uniprot Logical. If TRUE, UniProt protein IDs will be retrieved and included
#'        in the output. Default is TRUE.
#'
#' @return A data frame containing the following columns:
#'         - `UNIPROT`: UniProt protein IDs (if `get_uniprot` is TRUE).
#'         - `ENSEMBL_ID`: Ensembl gene IDs (if `get_ensembl` is TRUE).
#'         - `CHEBI_ID`: ChEBI IDs from the metabolomics data.
#'         The UniProt and Ensembl IDs are filtered to include only those present in
#'         `rowdata_target`. ChEBI IDs are repeated to match the lengths of the corresponding
#'         UniProt/Ensembl data.
#'
#' @details
#' The function performs the following steps:
#' 1. Extracts ChEBI IDs from the metabolomics data (`rowdata_metabo`).
#' 2. Filters out NA values from the ChEBI IDs.
#' 3. Queries UniProt for each ChEBI ID to retrieve associated protein and gene information.
#' 4. Filters the retrieved UniProt and Ensembl IDs to include only those present in
#'    the row names of `rowdata_target`.
#' 5. Constructs a data frame containing the relationships between ChEBI IDs, UniProt IDs,
#'    and Ensembl IDs.
#'
#'
#' @examples
#' # Example usage:
#' # relationships <- build_chebi_relatioships(rowdata_metabo, rowdata_target)
build_chebi_relationships <- function(rowdata_metabo, rowdata_target, get_ensembl = TRUE, get_uniprot = TRUE) {
  # Get the ChEBI IDs from the metabolomics data
  chebi_ids <- rownames(rowdata_metabo)

  # Filter out NA values
  chebi_ids <- chebi_ids[!is.na(chebi_ids)]

  # Query UniProt for the ChEBI IDs using the asynchronous function
  uniprot_data <- uniprot_chebi_query_async(chebi_id, rowdata_target, get_ensembl = get_ensembl, get_uniprot = get_uniprot)
  
  # Extract protein and gene IDs
  protein_ids <- unlist(lapply(uniprot_data, `[[`, "protein_ids"))
  gene_ids <- unlist(lapply(uniprot_data, `[[`, "gene_ids"))

  # Create a data frame with the results
  # This data frame includes UniProt IDs, Ensembl IDs, and ChEBI IDs
  # UniProt IDs and Ensembl IDs are filtered based on their presence in rowdata_target
  # ChEBI IDs are repeated to match the lengths of the corresponding UniProt/Ensembl data
  relationships_df <- data.frame(
    CHEBI_ID = rep(chebi_ids, ifelse(get_uniprot, lengths(uniprot_data$protein_ids), lengths(uniprot_data$gene_ids))),
  )

  if (get_uniprot) {
    relationships_df$UNIPROT <- protein_ids
  }

  if (get_ensembl) {
    relationships_df$ENSEMBL_ID <- gene_ids
  }

  return(relationships_df)
}

#' Build UniProt to Ensembl Mapping and Add ChEBI-Ensembl Relationships
#'
#' This script contains two functions:
#' 1. `build_uniprot_to_ensembl`: Maps UniProt IDs to Ensembl IDs for a given organism
#'    and filters the results based on the presence of Ensembl IDs in a provided transcript dataset.
#' 2. `add_chebi_ensembl_relationships`: Merges ChEBI-to-UniProt relationships with
#'    UniProt-to-Ensembl relationships to establish ChEBI-to-Ensembl mappings.
#'
#' ## Functions:
#'
#' ### `build_uniprot_to_ensembl(rowdata_prot, rowdata_trans, organism)`
#'
#' - **Purpose**: Maps UniProt protein IDs to Ensembl transcript IDs for a specified organism.
#' - **Parameters**:
#'   - `rowdata_prot`: A data frame where row names are UniProt IDs.
#'   - `rowdata_trans`: A data frame where row names are Ensembl transcript IDs.
#'   - `organism`: A string specifying the organism ("Hs" for human, "Mm" for mouse).
#' - **Details**:
#'   - Uses the appropriate organism annotation database (`org.Hs.eg.db` or `org.Mm.eg.db`).
#'   - Maps UniProt IDs to Ensembl IDs using the `AnnotationDbi::select` function.
#'   - Warns if any UniProt IDs cannot be mapped to Ensembl IDs.
#'   - Stops execution if no mappings are found.
#'   - Filters the results to include only Ensembl IDs present in `rowdata_trans`.
#' - **Returns**: A filtered data frame containing UniProt-to-Ensembl mappings.
#'
#' ### `add_chebi_ensembl_relationships(relationships_chebi_to_uniprot, relationships_uniprot_to_ensembl)`
#'
#' - **Purpose**: Combines ChEBI-to-UniProt relationships with UniProt-to-Ensembl mappings.
#' - **Parameters**:
#'   - `relationships_chebi_to_uniprot`: A data frame containing ChEBI-to-UniProt relationships.
#'   - `relationships_uniprot_to_ensembl`: A data frame containing UniProt-to-Ensembl mappings.
#' - **Details**:
#'   - Performs an outer join on the `UNIPROT` column to merge the two datasets.
#'   - Ensures that all relationships are preserved, even if some entries are missing in one dataset.
#' - **Returns**: A merged data frame containing ChEBI-to-Ensembl relationships.
build_uniprot_to_ensembl <- function(rowdata_prot, rowdata_trans, organism) {
  
  # Create a new column in the data frame based on the keys_list
  # Choose the appropriate annotation database based on the organism

  # Select the appropriate organism database
  if (organism == "Hs") {
    orgdb <- org.Hs.eg.db
  } else if (organism == "Mm") {
    orgdb <- org.Mm.eg.db
  } else {
    stop("Unsupported organism. Please use 'Hs' (human) or 'Mm' (mouse).")
  }

  # Define the keys list for mapping
  keys_list <- rownames(rowdata_prot)

  new_column <- AnnotationDbi::select(
    orgdb,
    keys = keys_list,
    columns = "ENSEMBL",
    keytype = "UNIPROT" # Assuming the source_type is UNIPROT
  )

  # Check for missing values and warn if any entries were not found
  num_missing <- sum(is.na(new_column$ENSEMBL))
  if (num_missing > 0) {
    warning(sprintf(
      "Mapping from %s to %s: %d/%d entries not found (returned NA).",
      "UNIPROT", "ENSEMBL", num_missing, length(keys_list)
    ))
  }

  # If all entries are missing, stop the function with an error
  if (num_missing == length(keys_list)) {
    stop(sprintf("Mapping from %s to %s failed: No entries could be mapped. Are you sure your UNIPROTS are valid?", "UNIPROT", "ENSEMBL"))
  }

  # Filter the rows to keep only those with an Ensembl ID present in rownames(rowdata_trans)
  return(new_column[new_column$ENSEMBL %in% rownames(rowdata_trans), ])
}

# NOT USED
add_chebi_ensembl_relationships <- function(relationships_chebi_to_uniprot, relationships_uniprot_to_ensembl) {
  # Perform an outer join on the UNIPROT column
  merged_relationships <- merge(
    relationships_chebi_to_uniprot,
    relationships_uniprot_to_ensembl,
    by = "UNIPROT",
    all = TRUE
  )

  return(merged_relationships)
}

chebi_relationships <- function(se_metabo, se_uniprot, organism, get_ensembl = FALSE, get_uniprot = FALSE) {
  # Throw an error if both get_ensembl and get_uniprot are FALSE
  if (!get_ensembl && !get_uniprot) stop("Both get_ensembl and get_uniprot cannot be FALSE. At least one must be TRUE.")

  # Build ChEBI relationships based on the provided SummarizedExperiment objects
  build_chebi_relationships(rowData(se_metabo), rowData(se_uniprot), get_ensembl = get_ensembl, get_uniprot = get_uniprot, organism)
}

uniprot_relationships <- function(se_uniprot, se_trans, organism) {
  build_uniprot_to_ensembl(rowData(se_uniprot), rowData(se_trans), organism)
}
