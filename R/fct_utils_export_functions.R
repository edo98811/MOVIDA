#' Export Code to File
#'
#' Function to be used in MOVIDA app to export code snippets to a file.
#' 
#' This function exports a given code object to a file, using a specified ID string.
#' It processes the code and ID, generating a formatted string that includes variable
#' definitions and the provided code. This is useful for saving reproducible code snippets
#' or analysis steps for later use or sharing.
export_code <- function(code, id) {

  # Convert the code object to a character string
  code_text <- paste(deparse(code)[-c(1, 2, length(deparse(code)))], collapse = "\n")

  # Remove "export_data = export_data" and the preceding comma from the code text
  code_text <- gsub(",\\s*export_data\\s*=\\s*export_data", "", code_text)

  # Split the ID string into individual assignments
  assignments <- strsplit(id, ";")[[1]]
  
  # Define the initial variable definitions
  var_def <- "(MOVIDA)\nmovida_data <- MOVIDA::movida_data(movida_list)\n\n"
  
  # Process each assignment and generate variable definitions
  var_def <- paste(var_def, paste(sapply(assignments, function(assignment) {
    # Split each assignment into variable name and value
    var_def_str <- strsplit(assignment, "_")[[1]]
    # Create a string for the variable definition
    paste0(var_def_str[[1]], " <- '", var_def_str[[2]], "'")
  }), collapse = "\n"), "\n")
  
  var_def <- paste(var_def, "\n", sep = "")

  # Combine the variable definitions and the code into the final output
  writeLines(paste(var_def, code_text, sep = "\n"), "output.txt")
  return(paste(var_def, code_text, sep = "\n"))
}


# export_code(code, id)
