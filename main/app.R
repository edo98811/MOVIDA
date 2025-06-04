MovidaApp <- function(movida_list) {
  # Load required package
  library(shiny)
  
  # Source any scripts your app depends on
  source("main/data_model.R")
  source("main/ui.R")
  source("main/server.R")
  
  # Initialize the global data model
  # assign("movida_data", MovidaModel$new(movida_list), envir = .GlobalEnv)
  assign("movida_data", MovidaModel$new(movida_list), envir = .GlobalEnv)
  
  # Wrap server to inject movida_data as argument
  server_wrapper <- function(input, output, session) {
    server(input, output, session, movida_data = movida_data)
  }
  
  shinyApp(ui = ui, server = server_wrapper)
}
