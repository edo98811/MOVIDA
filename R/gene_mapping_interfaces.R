
#' General gene ID mapping interface
#'  
#' This function provides a unified interface for mapping gene IDs across different databases.
#' @param from Source ID type (e.g., "ENSEMBL", "SYMBOL", "UNIPROT")
#' @param to Target ID type (e.g., "KEGG", "ENTREZID")
#' @param genes Vector of gene IDs to map
#' @param organism Organism name ("human" or "mouse")
#' @param db Optional AnnotationDbi database object. If NULL, it will be selected based on the organism from annotationHub.
#' @param ... Additional arguments passed to the annotationHub mapIds function.
#' @return Named vector of mapped IDs, names are original IDs
#' @export
gene_based_interface <- function(from, to, genes, organism = "human", db = NULL, ...) {

  # Handle the conversion based on the to parameter if to is KEGG first convert to ENTREZID then to KEGG
  res <- switch(
    to,
    KEGG = {
      entrez <- annotation_dbi_interface(from, "ENTREZID", genes, organism = organism, db = db, ...)
      entrez_to_kegg(entrez, organism = organism)
    },
    {
      annotation_dbi_interface(from, to, genes, organism = organism, db = db, ...)
    }
  )

return (res)


}
  #  Available databses: https://www.uniprot.org/help/id_mapping
  #    [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
  #  [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"
  # [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MGI"
  # [16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"
  # [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UNIPROT"

#' Select the appropriate AnnotationDbi database based on organism
#'
#' @param organism Organism name ("human" or "mouse")
#' @return The corresponding AnnotationDbi database
#' @importFrom AnnotationHub AnnotationHub query
#' @export
select_database <- function(organism) {
  ah <- AnnotationHub()
  if (organism %in% c("human", "hsa")) {
    db <- query(ah, c("Homo sapiens", "OrgDb"))[[1]]
  } else if (organism %in% c("mouse", "mmu")) {
    db <- query(ah, c("Mus musculus", "OrgDb"))[[1]]
  } else {
    stop("Unsupported organism. Please use 'human' or 'mouse'.")
  }
  return(db)
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
  if (is.null(db)) db <- select_database(organism)
  ids <- mapIds(db,
    keys = keys,
    column = to,
    keytype = from,
    ...
  )
  return(ids)
}
