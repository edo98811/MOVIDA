

library(webchem)

# Function to get ChEBI IDs from a list of InChIKeys
get_chebi_id_from_inchikey <- function(inchikey_list) {

  # Use lapply to iterate over each InChIKey in the list
  chebi_ids <- lapply(inchikey_list, function(inchikey) {
    tryCatch({
      # Attempt to retrieve the ChEBI ID using the webchem package
      get_chebi(inchikey, from = "inchikey")
    }, error = function(e) {
      # Handle errors gracefully and provide a message
      message(paste("Error retrieving ChEBI ID for InChIKey:", inchikey, "-", e$message))
      return(NA) # Return NA if an error occurs
    })
  })

  # Return the list of ChEBI IDs
  return(chebi_ids)
}
