# Function to check if rownames(rowData(se)) are of type Ensembl
check_ensembl <- function(names_to_check) {
  invalid_names <- names_to_check[!grepl("^ENS[A-Z]*[0-9]+$", names_to_check)]
  if (length(invalid_names) > 0) {
    warning("Warning: The following rownames(rowData(se)) are not valid Ensembl IDs: ", paste(invalid_names, collapse = ", "))
    return(FALSE)
  }
  return(TRUE)
}

# Function to check if rownames(rowData(se)) are of type Ensembl
check_goterm <- function(names_to_check) {
  invalid_names <- names_to_check[!grepl("^GO:[A-Z]*[0-9]+$", names_to_check)]
  if (length(invalid_names) > 0) {
    warning("Warning: The following rownames(rowData(se)) are not valid go-terms: ", paste(invalid_names, collapse = ", "))
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

check_symbol <- function(names_to_check) {

  # To make sure that the gene symbols are in lower case (harmonize human and mouse data)
  names_to_check <- tolower(names_to_check)
  invalid_names <- names_to_check[!grepl("^[a-z0-9]+$", names_to_check)]
  
  if (length(invalid_names) > 0) {
    warning("Warning: The following rownames(rowData(se)) are not valid gene symbols: ", paste(invalid_names, collapse = ", "))
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

check_metadata <- function(movida_list) {

  # Check that colnames of se_metabo, se_prot, se_trans and rownames of metadata are equal
  cn_metabo <- sort(colnames(movida_list$se_metabo))
  cn_prot <- sort(colnames(movida_list$se_prot))
  cn_trans <- sort(colnames(movida_list$se_trans))
  rn_metadata <- sort(rownames(movida_list$metadata))

  all_equal <- identical(cn_metabo, cn_prot) &&
               identical(cn_metabo, cn_trans) &&
               identical(cn_metabo, rn_metadata)

  if (!all_equal) {
    stop("Error: colnames of se_metabo, se_prot, se_trans and rownames of metadata must be identical (give as input only if experiment is matched and has same samples).")
  }
}

# Function to check if movida_list contains the required elements
check_movida_list <- function(movida_list) {

  # Check if at least one se_ object is not NULL
  if (all(sapply(movida_list[c("se_prot", "se_trans", "se_metabo")], is.null))) {
    stop("Error movida_list: At least one se_ object must not be NULL.")
  }

  # Check if at least one results_ object is not NULL
  # if (all(sapply(movida_list[c("results_prot", "results_trans", "results_metabo")], is.null))) {
  #   stop("Error movida_list: At least one results_object must not be NULL.")
  # }

  # Check if se_ objects inherit from SummarizedExperiment
  se_objects <- movida_list[c("se_prot", "se_trans", "se_metabo")]
  invalid_objects <- names(se_objects)[!sapply(se_objects, function(se) is.null(se) || inherits(se, "SummarizedExperiment"))]
  if (length(invalid_objects) > 0) {
    stop("Error movida_list: The following se_ objects do not inherit from SummarizedExperiment: ", paste(invalid_objects, collapse = ", "))
  }

  # Check if se_ objects are valid
  check_se(se_objects)

  # Check if organism is valid
  if (!check_organism(movida_list$organism)) {
    stop("Error movida_list: Invalid organism.")
  }

  return(TRUE)
}

check_movida_list_dde <- function(movida_list) {
  # Check if at least one se_ object is not NULL
  if (all(sapply(movida_list[c("dde_prot", "dde_trans", "dde_metabo")], is.null))) {
    stop("Error movida_list: At least one dde_object must not be NULL.")
  }

  # Check the objects that are of expected type
  dde_objects <- movida_list[c("dde_prot", "dde_trans", "dde_metabo")]
  invalid_objects <- names(dde_objects)[!sapply(dde_objects, function(se) is.null(se) || inherits(se, "SummarizedExperiment"))]
  if (length(invalid_objects) > 0) {
    stop("Error movida_list: The following dde_objects are not valid: ", paste(invalid_objects, collapse = ", "))
  }

  invalid_objects <- names(dde_objects)[!sapply(dde_objects, function(se) is.null(se) || inherits(se, "DeeDeeExperiment"))]
  if (length(invalid_objects) > 0) {
    warning("Error movida_list: The following dde_objects are ummarized experiment, the app will work but for full functionality use deedee experiments: ", paste(invalid_objects, collapse = ", "))
  }

  # Check metadata
  for (dde in dde_objects) {
    if (!is.null(dde)) {
      if (!"group" %in% colnames(colData(dde))) {
        stop("Error: colData of the SummarizedExperiment object must contain a 'group' column.")
      }
    }
  }

  # Check if groups overlap across se_objects
  group_lists <- lapply(dde_objects, function(dde) {
    if (!is.null(dde)) {
      return(unique(colData(dde)$group))
    } else {
      return(NULL)
    }
  })

  # Remove NULLs and check for overlap
  group_lists <- Filter(Negate(is.null), group_lists)
  if (length(group_lists) > 1 && length(Reduce(intersect, group_lists)) == 0) {
    warning("Warning: Groups do not overlap at all across se_objects.")
  }

  # Check if rownames(rowData(se)) are valid for each se_ object
  if (!is.null(dde_objects$dde_metabo) && !suppressWarnings(check_chebi(rownames(rowData(dde_objects$dde_metabo)))) && !suppressWarnings(check_inchi(rownames(rowData(se_objects$se_metabo))))) {
    stop("Error: Invalid ChEBI or InChI IDs in rownames(rowData(dde_metabo)).")
  }
  if (!is.null(dde_objects$dde_prot) && !suppressWarnings(check_uniprot(rownames(rowData(dde_objects$dde_prot))))) {
    stop("Error: Invalid UniProt IDs in rownames(rowData(dde_prot)).")
  }
  if (!is.null(dde_objects$dde_trans) && !suppressWarnings(check_ensembl(rownames(rowData(dde_objects$dde_trans)))) && !suppressWarnings(check_symbol(rownames(rowData(se_objects$se_trans))))) {
    stop("Error: Invalid Ensembl IDs or Gene Symbols in rownames(rowData(dde_trans)).")
  }

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

  # Check if rownames(rowData(se)) are valid for each se_ object
  if (!is.null(se_objects$se_metabo) && !suppressWarnings(check_chebi(rownames(rowData(se_objects$se_metabo)))) && !suppressWarnings(check_inchi(rownames(rowData(se_objects$se_metabo))))) {
    stop("Error: Invalid ChEBI or InChI IDs in rownames(rowData(se_metabo)).")
  }
  if (!is.null(se_objects$se_prot) && !suppressWarnings(check_uniprot(rownames(rowData(se_objects$se_prot))))) {
    stop("Error: Invalid UniProt IDs in rownames(rowData(se_prot)).")
  }
  if (!is.null(se_objects$se_trans) && !suppressWarnings(check_ensembl(rownames(rowData(se_objects$se_trans))))) {
    stop("Error: Invalid Ensembl IDs or Gene Symbols in rownames(rowData(se_trans)).")
  }

  return(TRUE)
}

# Function to validate contrasts
check_contrast <- function(contrast, groups) {
  split <- strsplit(contrast, "_vs_")[[1]]
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
