#' Convert KEGG pathway ID to readable pathway name
#'
#' @param pathway_id A KEGG pathway ID (e.g., "hsa04110", "mmu04110", "04110")
#' @param organism Optional organism code (e.g., "hsa", "mmu"). If not provided,
#'   it will be extracted from the pathway_id if present.
#' @return A character string with the readable pathway name
get_pathway_name <- function(id) {
  tryCatch(
    {
      res <- KEGGREST::keggGet(id)
      res[[1]]$NAME
    },
    error = function(e) {
      warning(paste("Could not retrieve pathway name for ID:", id))
      return("")
    }
  )
}

#' Download KEGG KGML file for a given pathway ID, with caching.
#'
#' @param pathway_id KEGG pathway ID (e.g., "hsa04110" or "04110").
#'
#' @return Content of the respose as text, or NULL if download failed.
get_kgml <- function(url) {
  # Download
  resp <- httr::GET(url)

  # Check if the request was successful (status 200)
  if (httr::http_error(resp)) {
    warning(
      "Failed to download KGML from URL: ", url,
      " (HTTP status ", httr::status_code(resp), ")"
    )
    return(NULL)
  }

  content <- httr::content(resp, as = "text", encoding = "UTF-8")
  return(content)
}

#' Get KEGG compounds with caching.
#' @param bfc A BiocFileCache object for caching.
#' @return A named character vector of KEGG compounds.
#' @importFrom BiocFileCache bfcrid bfcadd bfcpath
#' @importFrom KEGGREST keggList
get_compounds <- function(bfc) {

  # Check if cache exists
  cache_key <- "kegg_compound_cache"
  rid <- bfcrid(bfc, cache_key)

  # If cached, load from cache
  if (length(rid) > 0) {
    message("Loading KEGG compounds from cache...")
    compounds <- readRDS(bfcpath(bfc, rid))
    return(compounds)
  }

  # Otherwise download from KEGG
  message("Downloading KEGG compounds...")
  kegg_compounds <- KEGGREST::keggList("compound")

  # Save to cache
  rid <- bfcadd(bfc, cache_key, file = tempfile(fileext = ".rds"))
  saveRDS(kegg_compounds, bfcpath(bfc, rid))

  return(kegg_compounds)
}
