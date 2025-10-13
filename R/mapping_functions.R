#' Ensemble ID → KEGG ID
#' 
#' Map a vector of Ensembl gene IDs to KEGG gene IDs using UniProt as an intermediary
#' 
#' @param ensembl_ids Vector of Ensembl gene IDs
#' @param organism Organism name 
#' 
#' @return Named vector of KEGG gene IDs, names are original Ensembl IDs
#' 
#' @importFrom AnnotationDbi mapIds
#' 
#' @export
ensembl_to_kegg <- function(ensembl_ids, organism = "human") {
  uniprot_ids <- ensembl_to_uniprot(ensembl_ids, organism = organism)
  kegg_ids <- uniprot_to_kegg_batch(uniprot_ids, organism = organism)
  return(kegg_ids)
}

#' Gene Symbol → KEGG
#' 
#' Map a vector of gene symbols to KEGG gene IDs using UniProt as an intermediary
#' 
#' @param symbols Vector of gene symbols
#' @param organism Organism name
#' 
#' @return Named vector of KEGG gene IDs, names are original gene symbols
#' 
#' @importFrom AnnotationDbi mapIds
#' 
#' @export
symbol_to_kegg <- function(symbols, organism = "human") {
  uniprot_ids <- symbol_to_uniprot(symbols, organism = organism)
  kegg_ids <- uniprot_to_kegg_batch(uniprot_ids, organism = organism)
  return(kegg_ids)
}

uniprot_to_kegg_batch <- function(uniprot_ids, organism = "human", batch_size = 400) {
  organism <- to_organism_kegg(organism)

  all_results <- c()

  for (i in seq(1, length(uniprot_ids), by = batch_size)) {
    batch_ids <- uniprot_ids[i:min(i + batch_size - 1, length(uniprot_ids))]
    mapped_ids <- map_ids("UniProtKB_AC-ID", "KEGG", batch_ids)
    all_results <- c(all_results, mapped_ids)
  }

  # Ensure output order matches input
  all_results <- all_results[uniprot_ids]

  return(all_results)
}

#' Select the appropriate AnnotationDbi database based on organism
#'  
#'  @param organism Organism name ("human" or "mouse")
#' @return The corresponding AnnotationDbi database
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @return The corresponding AnnotationDbi database
#' @export
select_database <- function(organism) {
  if (organism %in% c("human", "hsa")) {
    return(org.Hs.eg.db)
  } else if (organism %in% c("mouse", "mmu")) {
    return(org.Mm.eg.db)
  } else {
    stop("Unsupported organism. Please use 'human' or 'mouse'.")
  }
}

#' UniProt ID → KEGG ID
#' Map a vector of UniProt IDs to KEGG gene IDs using UniProt ID mapping service
#' 
#' @param uniprot_ids Vector of UniProt IDs
#' @param organism Organism name
#' 
#' @return Named vector of KEGG gene IDs, names are original UniProt IDs
#'  
#' @export
uniprot_to_kegg <- function(...) {
  return(uniprot_to_kegg_batch(...))
}


#' Entrez ID → KEGG ID
#'
#' Map a vector of Entrez gene IDs to KEGG gene IDs using UniProt as an intermediary
#'
#' @param entrez_ids Vector of Entrez gene IDs
#' @param organism Organism name
#'
#' @return Named vector of KEGG gene IDs, names are original Entrez IDs
#'
#' @export
entrez_to_kegg <- function(entrez_ids, organism = "human") {
  # Convert organism name to KEGG code
  organism_code <- to_organism_kegg(organism)
  
  # Build KEGG-style IDs (e.g. "hsa:7157")
  kegg_ids <- paste0(organism_code, ":", entrez_ids)
  
  return(kegg_ids)
}

#' Map IDs using AnnotationDbi
#'
#' @param from Source ID type (e.g., "ENSEMBL", "SYMBOL", "   UNIPROT")
#' @param to Target ID type (e.g., "KEGG", "ENTREZID")
#' @param keys Vector of IDs to map
#' @param organism Organism name
#' @return Named vector of mapped IDs, names are original IDs
#' @importFrom AnnotationDbi mapIds
#' @export
annotation_dbi_interface <- function(from, to, keys, organism = "human", db = NULL, ...) {
  if(is.null(db)) db <- select_database(organism)
  ids <- mapIds(db,
    keys = keys,
    column = to,
    keytype = from,
    ...
  )
  return(ids)
}