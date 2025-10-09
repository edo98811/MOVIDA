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
