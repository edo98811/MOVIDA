  library(httr)
library(jsonlite)


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

build_uniprot_to_ensembl <- function(rowdata_prot, rowdata_trans, orgdb) {

  uniprot_ids <- rownames(rowdata_prot)

  # Create a new column in the data frame based on the keys_list
  new_column <- AnnotationDbi::select(
        orgdb,
        keys = keys_list,
        columns = "ENSEMBL",
        keytype = "UNIPROT"  # Assuming the source_type is UNIPROT
      )

  # Check for missing values and warn if any entries were not found
  num_missing <- sum(is.na(new_column))
  if (num_missing > 0) {
    warning(sprintf("Mapping from %s to %s: %d/%d entries not found (returned NA).",
                    source_type, destination_type, num_missing, length(keys_list)))
  }

  # If all entries are missing, stop the function with an error
  if (num_missing == length(keys_list)) {
    stop(sprintf("Mapping from %s to %s failed: No entries could be mapped. Are you sure your UNIPROTS are valid?")
  }

  # Filter the rows to keep only those with an Ensembl ID present in rownames(rowdata_trans)
  filtered_column <- new_column[new_column$ENSEMBL %in% rownames(rowdata_trans), ]

  return(new_column)
}

get_chebi_from_uniprot <- function(uniprot_id) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, ".json")
  res <- GET(url)
  stop_for_status(res)
  data <- fromJSON(content(res, "text"))
  
  chebi_ids <- c()

  # Look into ligand section for ChEBI references
  if(!is.null(data$comments)) {
    for(comment in data$comments) {
      if(comment$type == "cofactor" || comment$type == "binding") {
        for(ligand in comment$ligand) {
          if(!is.null(ligand$id) && grepl("^CHEBI:", ligand$id)) {
            chebi_ids <- c(chebi_ids, ligand$id)
          }
        }
      }
    }
  }
  unique(CHEBI_ID)
}

get_uniprot_from_chebi <- function(chebi_id, reviewed = TRUE, limit = 500) {
  base_url <- "https://rest.uniprot.org/uniprotkb/search"

  query <- paste0("ligand.id:", chebi_id)
  if (reviewed) {
    query <- paste0(query, " AND reviewed:true")
  }

  res <- GET(
    url = base_url,
    query = list(
      query = query,
      format = "json",
      size = limit,
      fields = "accession,protein_name,gene_names"
    )
  )

  stop_for_status(res)
  results <- fromJSON(content(res, "text"))
  
  if (length(results$results) == 0) {
    return(data.frame())
  }

  data.frame(
    UNIPROT_ID = sapply(results$results, function(x) x$primaryAccession),
    protein_name = sapply(results$results, function(x) x$proteinDescription$recommendedName$fullName$value),
    SYMBOL = sapply(results$results, function(x) paste(x$genes[[1]]$geneName$value, collapse = ","))
  )
}

chebi_to_uniprot <- function(chebi_ids, reviewed = TRUE, limit = 500) {
  if (length(chebi_ids) == 0) {
    return(data.frame())
  }
  
  results <- lapply(chebi_ids, get_uniprot_for_chebi, reviewed = reviewed, limit = limit)
  results_df <- do.call(rbind, results)
  
  if (nrow(results_df) == 0) {
    return(data.frame())
  }
  
  results_df <- unique(results_df)
  rownames(results_df) <- NULL
  return(results_df)
}

# The point is that in this case, I query the uniprot database for the uniprot that are associated with a secific ligand (chebi_id)
chebi_to_uniprot <- function(chebi_id) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/search?query=ligand-id:", chebi_id,
                "&fields=accession&format=json")
  res <- GET(url)
  data <- fromJSON(content(res, "text", encoding = "UTF-8"))
  if (length(data$results) == 0) return(NULL)
  sapply(data$results, function(x) x$primaryAccession)
}


uniprot_chebi_query <- function(chebi_id, get_ensembl = FALSE, get_uniprot = TRUE) {
  base_url <- "https://rest.uniprot.org/uniprotkb/search"

  params <- list(
    query = "29105 AND reviewed:true",
    fields = "accession,protein_name,cc_function,ft_binding",
    sort = "accession desc",
  )

  req <- request(base_url)
  req |> req_headers(
    accept = "application/json"
  )
  req |> req_url_query(!!!params)
  resp <- req_perform(req)

  if (resp_status(resp) != 200) {
    stop(sprintf("Error %d: %s", resp_status(resp), resp_body_string(resp)))

  data <- fromJSON(content(resp, "text", encoding = "UTF-8"))

  if (get_uniprot) {
    protein_ids <- sapply(data$results, function(x) x$primaryAccession)
  } else {
     protein_ids <- NULL
  }

  if (get_ensembl) {
    gene_ids <- sapply(data$results, function(x) {
      if (!is.null(x$uniProtKBCrossReferences)) {
        ensembl_ref <- Filter(function(ref) ref$database == "Ensembl", x$uniProtKBCrossReferences)
        if (length(ensembl_ref) > 0 && !is.null(ensembl_ref[[1]]$properties)) {
          gene_property <- Filter(function(prop) prop$key == "GeneId", ensembl_ref[[1]]$properties)
          if (length(gene_property) > 0) {
            return(gene_property[[1]]$value)
          }
        }
      }
      return(NA)
    })
  } else {
    gene_ids <- NULL
  }

  return(list(protein_ids = protein_ids, gene_ids = gene_ids))
}

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

