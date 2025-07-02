#' @name MovidaApp
#' @title Launch the MOVIDA Shiny Application
#'
#' @description
#' This function initializes and launches the MOVIDA Shiny application using the provided list of MOVIDA data.
#' It loads required scripts, initializes the global data model, and starts the Shiny app.
#'
#' @param movida_list A list containing the data required to initialize the MOVIDA data model.
#'
#' @details
#' The function sources UI, server, and data model scripts, initializes the global \code{movida_data} object,
#' loads relationships from the \code{relationships/} directory, and launches the Shiny application.
#'
#' @return The function starts the Shiny application and does not return any value.
#'
#' @import shiny
#' @import bslib
#' @import DT
#' @export
MovidaApp <- function(movida_list) {

  # Initialize the global data model
  # assign("movida_data", MovidaModel$new(movida_list), envir = .GlobalEnv)
  assign("movida_data", MovidaModel$new(movida_list), envir = .GlobalEnv)
  movida_data$load_relationships("relationships/")
  
  # Wrap server to inject movida_data as argument
  server_wrapper <- function(input, output, session) {
    server(input, output, session, movida_data = movida_data)
  }
  
  shinyApp(ui = ui, server = server_wrapper)
  # shiny::runApp(list(ui = ui, server = server_wrapper), display.mode="showcase")
}
