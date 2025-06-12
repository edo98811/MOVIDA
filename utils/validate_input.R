# Function to check if rownames(rowData(se)) are of type Ensembl
check_ensembl <- function(names_to_check) {
  invalid_names <- names_to_check[!grepl("^ENS[A-Z]*[0-9]+$", names_to_check)]
  if (length(invalid_names) > 0) {
    warning("Warning: The following rownames(rowData(se)) are not valid Ensembl IDs: ", paste(invalid_names, collapse = ", "))
    return(FALSE)
  }
  return(TRUE)
}

# Function to check if rownames(rowData(se)) are of type ChEBI
check_chebi <- function(names_to_check) {
  invalid_names <- names_to_check[!grepl("^CHEBI:[0-9]+$", names_to_check)]
  if (length(invalid_names) > 0) {
    warning("Warning: The following rownames(rowData(se)) are not valid ChEBI IDs: ", paste(invalid_names, collapse = ", "))
    return(FALSE)
  }
  return(TRUE)
}

# Function to check if rownames(rowData(se)) are of type UniProt
check_uniprot <- function(names_to_check) {
  invalid_names <- names_to_check[!grepl("^[A-Z0-9]{6,10}$", names_to_check)]
  if (length(invalid_names) > 0) {
    warning("Warning: The following rownames(rowData(se)) are not valid UniProt IDs: ", paste(invalid_names, collapse = ", "))
    return(FALSE)
  }
  return(TRUE)
}

# Function to check if rownames(rowData(se)) are of type InChI
check_inchi <- function(names_to_check) {
  invalid_names <- names_to_check[!grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$", names_to_check)]
  if (length(invalid_names) > 0) {
    warning("Warning: The following rownames(rowData(se)) are not valid InChI identifiers: ", paste(invalid_names, collapse = ", "))
    return(FALSE)
  }
  return(TRUE)
}

# Function to check if the organism is either 'hs' or 'Mm'
check_organism <- function(organism) {
  if (!organism %in% c("Hs", "Mm")) {
    warning("Warning: organism must be either 'Hs' (Homo sapiens) or 'Mm' (Mus musculus).")
    return(FALSE)
  }
  return(TRUE)
}

# Function to check if movida_list contains the required elements
check_movida_list <- function(movida_list) {

  required_elements <- c(
    "results_prot", "results_trans", "results_metabo",
    "se_prot", "se_trans", "se_metabo", "organism"
  )
  # Check if at least one se_ object is not NULL
  if (all(sapply(movida_list[c("se_prot", "se_trans", "se_metabo")], is.null))) {
    stop("Error: At least one se_ object must not be NULL.")
  }

  # Check if at least one results_ object is not NULL
  if (all(sapply(movida_list[c("results_prot", "results_trans", "results_metabo")], is.null))) {
    stop("Error: At least one results_ object must not be NULL.")
  }

  # Check if se_ objects inherit from SummarizedExperiment
  se_objects <- movida_list[c("se_prot", "se_trans", "se_metabo")]
  if (!all(sapply(se_objects, function(se) is.null(se) || inherits(se, "SummarizedExperiment")))) {
    stop("Error: All se_ objects must either be NULL or inherit from SummarizedExperiment.")
  }

  # Check if colData of all se_ objects is valid
  if (!check_se(se_objects)) {
    stop("Error: colData of the SummarizedExperiment objects is not valid.")
  }

  # Check if organism is valid
  if (!check_organism(movida_list$organism)) {
    stop("Error: Invalid organism.")
  }

  return(TRUE)
}

# Function to check if colData of se objects contains a 'group' column
check_se <- function(se_objects) {
  for (se in se_objects) {
    if (!is.null(se)) {
      if (!"group" %in% colnames(colData(se))) {
        stop("Error: colData of the SummarizedExperiment object must contain a 'group' column.")
      }
    }
  }
  # Check if groups overlap across se_ objects
  group_lists <- lapply(se_objects, function(se) {
    if (!is.null(se)) {
      return(unique(colData(se)$group))
    } else {
      return(NULL)
    }
  })

  # Remove NULLs and check for overlap
  group_lists <- Filter(Negate(is.null), group_lists)
  if (length(group_lists) > 1 && length(Reduce(intersect, group_lists)) == 0) {
    warning("Warning: Groups do not overlap at all across se_ objects.")
  }
  return(TRUE)
}

# Function to validate contrasts
check_contrast <- function(contrast, groups) {
  split <- strsplit(contrast, "vs")[[1]]
  if (!all(split %in% groups)) {
    warning("Warning: Contrast contains invalid group names: ", contrast)
    return(FALSE)
  } else {
    return(TRUE)
  }
}

validate_contrasts <- function(contrasts, groups) {
  valid_contrasts <- contrasts[sapply(contrasts, function(contrast) check_contrast(contrast, groups))]
  return(valid_contrasts)
}
