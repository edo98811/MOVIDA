# Function to check if names_to_check are of type Ensembl
check_ensembl <- function(names_to_check) {
  invalid_names <- names_to_check[!grepl("^ENS[A-Z]*[0-9]+$", names_to_check)]
  if (length(invalid_names) > 0) {
    warning("check_ensembl: The following names_to_check are not valid Ensembl IDs: ", paste(invalid_names, collapse = ", "))
    return(FALSE)
  }
  return(TRUE)
}

# Function to check if names_to_check are of type Ensembl
check_goterm <- function(names_to_check) {
  invalid_names <- names_to_check[!grepl("^GO:[A-Z]*[0-9]+$", names_to_check)]
  if (length(invalid_names) > 0) {
    warning("check_goterm: The following names_to_check are not valid go-terms: ", paste(invalid_names, collapse = ", "))
    return(FALSE)
  }
  return(TRUE)
}

# Function to check if names_to_check are of type ChEBI
check_chebi <- function(names_to_check) {
  invalid_names <- names_to_check[!grepl("^CHEBI:[0-9]+$", names_to_check)]
  if (length(invalid_names) > 0) {
    warning("check_chebi: The following names_to_check are not valid ChEBI IDs: ", paste(invalid_names, collapse = ", "))
    return(FALSE)
  }
  return(TRUE)
}

# Function to check if names_to_check are of type UniProt
check_uniprot <- function(names_to_check) {
  invalid_names <- names_to_check[!grepl("^[A-Z0-9]{6,10}$", names_to_check)]
  if (length(invalid_names) > 0) {
    warning("check_uniprot: The following names_to_check are not valid UniProt IDs: ", paste(invalid_names, collapse = ", "))
    return(FALSE)
  }
  return(TRUE)
}

# Function to check if names_to_check are of type KEGG
check_kegg <- function(names_to_check) {
  invalid_names <- names_to_check[!grepl("^[A-Za-z]{2,4}:[0-9]+$", names_to_check)]
  if (length(invalid_names) > 0) {
    warning("check_kegg: The following names_to_check are not valid KEGG IDs: ", paste(invalid_names, collapse = ", "))
    return(FALSE)
  }
  return(TRUE)
}

# Function to check if names_to_check are of type InChI
check_inchi <- function(names_to_check) {
  invalid_names <- names_to_check[!grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$", names_to_check)]
  if (length(invalid_names) > 0) {
    warning("check_inchi: The following names_to_check are not valid InChI identifiers: ", paste(invalid_names, collapse = ", "))
    return(FALSE)
  }
  return(TRUE)
}

check_symbol <- function(names_to_check) {
  # To make sure that the gene symbols are in lower case (harmonize human and mouse data)
  names_to_check <- tolower(names_to_check)
  invalid_names <- names_to_check[!grepl("^[a-z0-9]+$", names_to_check)]

  if (length(invalid_names) > 0) {
    warning("check_symbol: The following names_to_check are not valid gene symbols: ", paste(invalid_names, collapse = ", "))
    return(FALSE)
  }
  return(TRUE)
}

# Function to check if names_to_check are of a valid feature type (supported types)
check_is_valid_feature <- function(names_to_check) {
  if (
      suppressWarnings(check_ensembl(names_to_check))
      || suppressWarnings(check_symbol(names_to_check))
      || suppressWarnings(check_uniprot(names_to_check))
      || suppressWarnings(check_chebi(names_to_check))
      || suppressWarnings(check_inchi(names_to_check))
      || suppressWarnings(check_kegg(names_to_check))
  ) {
    return(TRUE)
  } else {
    warning("check_is_valid_feature: Invalid source type.")
    return(FALSE)
  }
}

# Function to check if the organism is either 'hs' or 'Mm'
check_organism <- function(organism) {
  if (!organism %in% c("Hs", "Mm")) {
    warning("check_organism: organism must be either 'Hs' (Homo sapiens) or 'Mm' (Mus musculus).")
    return(FALSE)
  }
  return(TRUE)
}

# If common metadata is provided this checks that it correctly overlaps
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
    stop("check_metadata: colnames of se_metabo, se_prot, se_trans and rownames of metadata must be identical (give as input only if experiment is matched and has same samples).")
  }
}

# Checks that the inputs are valid objects
check_movida_list_dde <- function(movida_list) {
  # Check if at least one se_ object is not NULL
  if (all(sapply(movida_list[c("dde_prot", "dde_trans", "dde_metabo")], is.null))) {
    stop("check_movida_list_dde: At least one dde_object must not be NULL.")
  }

  # Check the objects that are of expected type
  dde_objects <- movida_list[c("dde_prot", "dde_trans", "dde_metabo")]
  invalid_objects <- names(dde_objects)[!sapply(dde_objects, function(se) is.null(se) || inherits(se, "SummarizedExperiment"))]
  if (length(invalid_objects) > 0) {
    stop("check_movida_list_dde: The following dde_objects are not valid: ", paste(invalid_objects, collapse = ", "))
  }

  invalid_objects <- names(dde_objects)[!sapply(dde_objects, function(se) is.null(se) || inherits(se, "DeeDeeExperiment"))]
  if (length(invalid_objects) > 0) {
    warning("check_movida_list_dde: The following dde_objects are summarized experiment, the app will work but for full functionality use deedee experiments: ", paste(invalid_objects, collapse = ", "))
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
    warning("check_movida_list_dde: Groups do not overlap at all across dde_objects")
  }
}

# Function to validate contrasts
check_contrast <- function(contrast, groups) {
  split <- strsplit(contrast, "_vs_")[[1]]
  if (!all(split %in% groups)) {
    warning("check_contrast: Contrast contains invalid group names: ", contrast)
    return(FALSE)
  } else {
    return(TRUE)
  }
}

validate_contrasts <- function(contrasts, groups) {
  valid_contrasts <- contrasts[sapply(contrasts, function(contrast) check_contrast(contrast, groups))]
  return(valid_contrasts)
}
