# Function to nicely format the size of an object in MB
format_object_size <- function(obj, unit = "MB") {
  size_in_bytes <- as.integer(object.size(obj))
  
  size_in_unit <- switch(
    unit,
    "B" = size_in_bytes,
    "KB" = size_in_bytes / 1000,
    "MB" = size_in_bytes / (1000 * 1000),
    "GB" = size_in_bytes / (1000 * 1000 * 1000),
    stop("Invalid unit. Choose from 'B', 'KB', 'MB', or 'GB'.")
  )
  
  formatted_size <- sprintf("%.2f %s", size_in_unit, unit)
  return(formatted_size)
}