  # Import relevant libraries
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)


build_chebi_relatioships <- function(rowdata_metabo, rowdata_second, get_ensembl = TRUE, get_uniprot = TRUE) {
  # Get the ChEBI IDs from the metabolomics data
  chebi_ids <- rowdata(rowData(se_metabo))

  # Filter out NA values
  chebi_ids <- chebi_ids[!is.na(chebi_ids)]

  # Query UniProt for the ChEBI IDs
  uniprot_data <- lapply(chebi_ids, function(chebi_id) {
    query_result <- uniprot_chebi_query(chebi_id, get_ensembl = get_ensembl, get_uniprot = get_uniprot)

    if (get_uniprot && !is.null(query_result$protein_ids)) {
      # Filter the UniProt IDs to include only those present in the rownames of rowdata_second
      query_result$protein_ids <- query_result$protein_ids[query_result$protein_ids %in% rownames(rowdata_second)]
    }

    if (get_ensembl && !is.null(query_result$gene_ids)) {
      # Filter the Ensembl IDs to include only those present in the rownames of rowdata_second
      query_result$gene_ids <- query_result$gene_ids[query_result$gene_ids %in% rownames(rowdata_second)]
    }

    query_result
  })

  # Extract protein and gene IDs
  protein_ids <- unlist(lapply(uniprot_data, `[[`, "protein_ids"))
  gene_ids <- unlist(lapply(uniprot_data, `[[`, "gene_ids"))

  # Create a data frame with the results
  # This data frame includes UniProt IDs, Ensembl IDs, and ChEBI IDs
  # UniProt IDs and Ensembl IDs are filtered based on their presence in rowdata_second
  # ChEBI IDs are repeated to match the lengths of the corresponding UniProt/Ensembl data
  result_df <- data.frame(
    UNIPROT_ID = if (get_uniprot) protein_ids else NA,
    ENSEMBL_ID = if (get_ensembl) gene_ids else NA,
    CHEBI_ID = rep(chebi_ids, lengths(uniprot_data))
  )

  return(result_df)
}


build_uniprot_to_ensembl <- function(rowdata_prot, rowdata_trans, organism) {
  
  uniprot_ids <- rownames(rowdata_prot)
  
  # Create a new column in the data frame based on the keys_list
  # Choose the appropriate annotation database based on the organism

  # Select the appropriate organism database
  if (organism == "Hs") {
    orgdb <- org.Hs.eg.db
  } else if (organism == "Mm") {
    orgdb <- org.Mm.eg.db
  }

  # Define the keys list for mapping
  keys_list <- uniprot_ids

  new_column <- AnnotationDbi::select(
    orgdb,
    keys = keys_list,
    columns = "ENSEMBL",
    keytype = "UNIPROT" # Assuming the source_type is UNIPROT
  )

  # Check for missing values and warn if any entries were not found
  num_missing <- sum(is.na(new_column))
  if (num_missing > 0) {
    warning(sprintf(
      "Mapping from %s to %s: %d/%d entries not found (returned NA).",
      "UNIPROT",  "ENSEMBL", num_missing, length(keys_list)
    ))
  }

  # If all entries are missing, stop the function with an error
  if (num_missing == length(keys_list)) {
    stop(sprintf("Mapping from %s to %s failed: No entries could be mapped. Are you sure your UNIPROTS are valid?"))
  }

  # Filter the rows to keep only those with an Ensembl ID present in rownames(rowdata_trans)
  filtered_column <- new_column[new_column$ENSEMBL %in% rownames(rowdata_trans), ]

  return(new_column)
}
