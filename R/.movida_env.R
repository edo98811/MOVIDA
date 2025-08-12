.movida_env <- new.env(parent = emptyenv())

get_movida_data <- function() {
  if (is.null(.movida_env$movida_data)) {
    stop("movida_data has not been initialized yet!")
  }
  .movida_env$movida_data
}


set_movida_data <- function(data) {
  .movida_env$movida_data <- data
}

initialize_movida_data <- function(movida_list) {
  movida_data <- MovidaModel$new(movida_list)
  movida_data$load_relationships("relationships/")
  set_movida_data(movida_data)
}