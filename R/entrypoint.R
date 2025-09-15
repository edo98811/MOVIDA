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
#' and launches the Shiny application.
#'
#' @return The function starts the Shiny application and does not return any value.
#'
#' @import shiny
#' @import bslib
#' @import DT
#' @export
MovidaApp <- function(movida_list, recording_test = FALSE) {

  movida_data <- MovidaModel$new(movida_list)
  movida_data$load_relationships("relationships/")

  # initialize_movida_data(movida_list)
  server_wrapper <- function(input, output, session) {
    server(input, output, session, movida_data = movida_data)
  }

  if (recording_test) {
    shinytest2::record_test(shinyApp(ui = ui, server = server_wrapper))
  } else {
    shinyApp(ui = ui, server = server_wrapper)
  }

}
