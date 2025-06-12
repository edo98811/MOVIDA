# Import relevant libraries
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(httr2)
library(jsonlite)
library(future)
library(furrr)
library(progressr)

# List of functions in this file:
# 1. uniprot_inchi_query_async: Asynchronously queries UniProt for multiple inchi IDs and filters results based on target data.
# 2. uniprot_inchi_query: Queries UniProt for a single inchi ID and retrieves associated protein and gene IDs.
# 3. build_inchi_relationships: Builds relationships between inchi IDs, UniProt IDs, and Ensembl IDs based on metabolomics and target data.
# 4. build_uniprot_to_ensembl: Maps UniProt IDs to Ensembl IDs for a given organism and filters results based on transcript data.
# 5. add_inchi_ensembl_relationships: NOT USED Merges inchi-to-UniProt relationships with UniProt-to-Ensembl mappings to establish inchi-to-Ensembl relationships.
# 6. inchi_relationships: Generates inchi relationships using SummarizedExperiment objects for metabolomics and target data.
# 7. uniprot_relationships: Builds UniProt-to-Ensembl relationships using SummarizedExperiment objects for protein and transcript data.

# Wrapper for async querying multiple inchi IDs
uniprot_inchi_query_async <- function(inchikeys, rowdata_target, get_ensembl = FALSE, get_uniprot = TRUE, sleep_time = 0.01, progressbar = NULL) {
  with_progress({
    p <- progressor(steps = length(inchikeys))
    results <- furrr::future_map(inchikeys, function(id) { # like purrr map but works in parallel: https://www.rdocumentation.org/packages/furrr/versions/0.3.1/topics/future_map
      p()
      # Sys.sleep(sleep_time) # Avoid hitting the API too hard (needed?)
      tryCatch(
        {
          # Query UniProt for each inchi ID (fucntion return list)
          query_result <- uniprot_inchi_query(id, get_ensembl = get_ensembl, get_uniprot = get_uniprot)

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
    },
    .options = furrr_options(seed = TRUE)
    ) # Ensures reproducibility
  })
  names(results) <- inchikeys
  return(results)
}

uniprot_inchi_query_sync <- function(inchikeys, rowdata_target, get_ensembl = FALSE, get_uniprot = TRUE, sleep_time = 0.01) {
  with_progress({
    p <- progressor(steps = length(inchikeys))
    results <- lapply(inchikeys, function(id) { # Sequential processing using lapply
      p()
      # Sys.sleep(sleep_time) # Avoid hitting the API too hard
      tryCatch(
        {
          # Query UniProt for each inchi ID (function returns list)
          query_result <- uniprot_inchi_query(id, get_ensembl = get_ensembl, get_uniprot = get_uniprot)

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
    })
  })
  names(results) <- inchikeys
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
# Function to query UniProt for inchi-related data
uniprot_inchi_query <- function(inchikey, get_ensembl = FALSE, get_uniprot = TRUE) {
  # Check if the provided inchi ID is valid
  if (!check_inchi(inchikey)) {
    stop("Not a valid inchi ID: ", inchikey)
  }

  # Base URL for the UniProt REST API
  base_url <- "https://rest.uniprot.org/uniprotkb/stream"

  # Define the query parameters for the API request
  # test query for postman: https://rest.uniprot.org/uniprotkb/stream?query=inchikey%3AKJTLQQUUPVSXIM-UHFFFAOYSA-N%20AND%20reviewed%3Atrue&fields=accession%2Cxref_ensembl&sort=accession%20desc
  params <- list(
    query = paste0("inchikey:", inchikey, " AND reviewed:true"),
    fields = paste(c(
      if (get_uniprot) "accession",
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
  data <- resp_body_json(resp)
  if (length(data$results) == 0) {
    warning("No matches found for the provided inchi ID:", inchikey)
    return(list(protein_ids = NULL, gene_ids = NULL))
    # stop("No matches found for the provided inchi ID during mapping to UniProt.")
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
    protein_ids <- unique(unlist(protein_ids)) # Convert to a character vector
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
    gene_ids <- unique(unlist(gene_ids[!is.na(gene_ids)])) # Convert to a character vector and remove NAs
  } else {
    gene_ids <- NULL
  }

  # Return a list containing the extracted protein and gene IDs
  return(list(protein_ids = protein_ids, gene_ids = gene_ids))
}

#' Build inchi Relationships
#'
#' This function generates a data frame of relationships between inchi IDs, UniProt IDs,
#' and Ensembl IDs based on metabolomics and target data. It queries UniProt to retrieve
#' protein and gene information associated with the provided inchi IDs and filters the
#' results based on their presence in the target data.
#'
#' @param rowdata_metabo A data frame or similar object containing metabolomics data.
#'        This should include inchi IDs.
#' @param rowdata_target A data frame or similar object containing target data.
#'        The row names should correspond to UniProt or Ensembl IDs.
#' @param get_ensembl Logical. If TRUE, Ensembl gene IDs will be retrieved and included
#'        in the output. Default is TRUE.
#' @param get_uniprot Logical. If TRUE, UniProt protein IDs will be retrieved and included
#'        in the output. Default is TRUE.
#'
#' @return A data frame containing the following columns:
#'         - `UNIPROT`: UniProt protein IDs (if `get_uniprot` is TRUE).
#'         - `ENSEMBL`: Ensembl gene IDs (if `get_ensembl` is TRUE).
#'         - `inchikey`: inchi IDs from the metabolomics data.
#'         The UniProt and Ensembl IDs are filtered to include only those present in
#'         `rowdata_target`. inchi IDs are repeated to match the lengths of the corresponding
#'         UniProt/Ensembl data.
#'
#' @details
#' The function performs the following steps:
#' 1. Extracts inchi IDs from the metabolomics data (`rowdata_metabo`).
#' 2. Filters out NA values from the inchi IDs.
#' 3. Queries UniProt for each inchi ID to retrieve associated protein and gene information.
#' 4. Filters the retrieved UniProt and Ensembl IDs to include only those present in
#'    the row names of `rowdata_target`.
#' 5. Constructs a data frame containing the relationships between inchi IDs, UniProt IDs,
#'    and Ensembl IDs.
#'
#'
#' @examples
#' # Example usage:
#' # relationships <- build_inchi_relatioships(rowdata_metabo, rowdata_target)
build_inchi_relationships <- function(rowdata_metabo, rowdata_target, get_ensembl = FALSE, get_uniprot = FALSE, progressbar = NULL) {

    # Throw an error if both get_ensembl and get_uniprot are FALSE
  if (!get_ensembl && !get_uniprot) stop("Both get_ensembl and get_uniprot cannot be FALSE. At least one must be TRUE.")

  # Check if get_uniprot or get_ensembl is TRUE and run check_type accordingly
  if (get_uniprot) {
    if (!check_uniprot(rownames(rowdata_target))) {
      stop("UniProt check failed: Invalid UniProt IDs in the target data.")
    }
  }

  if (get_ensembl) {
    if (!check_ensembl(rownames(rowdata_target))) {
      stop("Ensembl check failed: Invalid Ensembl IDs in the target data.")
    }
  }

  # Get the inchi IDs from the metabolomics data
  inchikeys <- rownames(rowdata_metabo)

  # Filter out NA values
  if (anyNA(inchikeys)) {
    warning("Some InChIKeys are NA and will be removed.")
    inchikeys <- inchikeys[!is.na(inchikeys)]
  }

  # Query UniProt for the InChIKeys
  uniprot_data <- uniprot_inchi_query_async(inchikeys, rowdata_target, get_ensembl = get_ensembl, get_uniprot = get_uniprot, progressbar = progressbar)

  # Extract IDs, handling missing lists safely
  protein_ids_list <- lapply(uniprot_data, function(x) x$protein_ids %||% character(0))
  gene_ids_list <- lapply(uniprot_data, function(x) x$gene_ids %||% character(0))

  # Determine how many times to repeat each InChIKey
  repeats <- if (get_uniprot) {
    lengths(protein_ids_list)
  } else {
    lengths(gene_ids_list)
  }

  # Create a data frame to store the relationships
  # Repeat each InChIKey based on the number of associated IDs (UniProt or Ensembl)
  relationships_df <- data.frame(
    INCHIKEY = rep(names(uniprot_data), repeats),
    stringsAsFactors = FALSE
  )

  # Append UniProt IDs to the data frame if requested
  if (get_uniprot) {
    relationships_df$UNIPROT <- unlist(protein_ids_list, use.names = FALSE)
  }

  # Append Ensembl IDs to the data frame if requested
  if (get_ensembl) {
    relationships_df$ENSEMBL <- unlist(gene_ids_list, use.names = FALSE)
  }

  return(relationships_df)
}

#' Build UniProt to Ensembl Mapping and Add inchi-Ensembl Relationships
#'
#' This script contains two functions:
#' 1. `build_uniprot_to_ensembl`: Maps UniProt IDs to Ensembl IDs for a given organism
#'    and filters the results based on the presence of Ensembl IDs in a provided transcript dataset.
#' 2. `add_inchi_ensembl_relationships`: Merges inchi-to-UniProt relationships with
#'    UniProt-to-Ensembl relationships to establish inchi-to-Ensembl mappings.
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
#' ### `add_inchi_ensembl_relationships(relationships_inchi_to_uniprot, relationships_uniprot_to_ensembl)`
#'
#' - **Purpose**: Combines inchi-to-UniProt relationships with UniProt-to-Ensembl mappings.
#' - **Parameters**:
#'   - `relationships_inchi_to_uniprot`: A data frame containing inchi-to-UniProt relationships.
#'   - `relationships_uniprot_to_ensembl`: A data frame containing UniProt-to-Ensembl mappings.
#' - **Details**:
#'   - Performs an outer join on the `UNIPROT` column to merge the two datasets.
#'   - Ensures that all relationships are preserved, even if some entries are missing in one dataset.
#' - **Returns**: A merged data frame containing inchi-to-Ensembl relationships.
build_uniprot_to_ensembl <- function(rowdata_prot, rowdata_trans, organism) {
  # Create a new column in the data frame based on the keys_list
  # Choose the appropriate annotation database based on the organism

  # Define the keys list for mapping
  keys_list <- rownames(rowdata_prot)

  orgdb <- switch(
    organism,
    "Hs" = org.Hs.eg.db,
    "Mm" = org.Mm.eg.db,
    stop("Unsupported organism. Please use 'Hs' for human or 'Mm' for mouse.")
  )
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
add_inchi_ensembl_relationships <- function(relationships_inchi_to_uniprot, relationships_uniprot_to_ensembl) {
  # Perform an outer join on the UNIPROT column
  merged_relationships <- merge(
    relationships_inchi_to_uniprot,
    relationships_uniprot_to_ensembl,
    by = "UNIPROT",
    all = TRUE
  )

  return(merged_relationships)
}

inchi_relationships <- function(se_metabo, se_target, organism, get_ensembl = FALSE, get_uniprot = FALSE) {

  # Throw an error if both get_ensembl and get_uniprot are FALSE
  if (!get_ensembl && !get_uniprot) stop("Both get_ensembl and get_uniprot cannot be FALSE. At least one must be TRUE.")
  # Check if get_uniprot or get_ensembl is TRUE and run check_type accordingly
  if (get_uniprot) {
    if (!check_uniprot(rownames(rowData(se_target)))) {
      stop("UniProt check failed: Invalid UniProt IDs in the target data.")
    }
  }

  if (get_ensembl) {
    if (!check_ensembl(rownames(rowData(se_target)))) {
      stop("Ensembl check failed: Invalid Ensembl IDs in the target data.")
    }
  }

  # Build inchi relationships based on the provided SummarizedExperiment objects
  build_inchi_relationships(rowData(se_metabo), rowData(se_target), get_ensembl = get_ensembl, get_uniprot = get_uniprot, organism)
}

uniprot_relationships <- function(se_uniprot, se_trans, organism) {
  build_uniprot_to_ensembl(rowData(se_uniprot), rowData(se_trans), organism)
}
