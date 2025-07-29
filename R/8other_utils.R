get_fea_list_wrapper <- function(dde, ...) {
  tryCatch({return(get_fea_list(dde, ...))},
  error = {return(NULL)}
  )
}

get_dea_list_wrapper <- function(dde, ...) {
  tryCatch({return(get_dea_list(dde, ...))},
  error = {return(NULL)}
  )
}

