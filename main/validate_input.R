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