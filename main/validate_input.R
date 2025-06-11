# Function to check if rownames(rowData(se)) are of type Ensembl
check_ensembl <- function(se) {
  if (!all(grepl("^ENS[A-Z]*[0-9]+$", rownames(rowData(se))))) {
    stop("Error: rownames(rowData(se)) are not all valid Ensembl IDs.")
  }
  return(TRUE)
}

# Function to check if rownames(rowData(se)) are of type ChEBI
check_chebi <- function(se) {
  if (!all(grepl("^CHEBI:[0-9]+$", rownames(rowData(se))))) {
    stop("Error: rownames(rowData(se)) are not all valid ChEBI IDs.")
  }
  return(TRUE)
}

# Function to check if rownames(rowData(se)) are of type UniProt
check_uniprot <- function(se) {
  if (!all(grepl("^[A-Z0-9]{6,10}$", rownames(rowData(se))))) {
    stop("Error: rownames(rowData(se)) are not all valid UniProt IDs.")
  }
  return(TRUE)
}

# Function to check if the organism is either 'hs' or 'Mm'
check_organism <- function(organism) {
  if (!organism %in% c("Hs", "Mm")) {
    stop("Error: organism must be either 'Hs' (Homo sapiens) or 'Mm' (Mus musculus).")
  }
  return(TRUE)
}

# Function to check if movida_list contains the required elements
check_movida_list <- function(movida_list) {
  required_elements <- c(
    "results_prot", "results_trans", "results_metabo",
    "se_prot", "se_trans", "se_metabo", "organism"
  )
  
  missing_elements <- setdiff(required_elements, names(movida_list))
  
  if (length(missing_elements) > 0) {
    stop(paste("Error: movida_list is missing the following elements:", 
               paste(missing_elements, collapse = ", ")))
  }
  
  return(TRUE)
}