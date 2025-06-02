library(shiny)
library(gargoyle)

# UI function for the module
mod_overview_ui <- function(id) {
  ns <- NS(id) # Namespace for the module
  tagList(
    h3("Overview Module Placeholder"),
    actionButton(ns("toggle_sidebar_btn"), "Toggle Sidebar")
  )
}

# Server function for the module
mod_overview_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$toggle_sidebar_btn, {
      trigger("toggle_sidebar")
    })
  })
}
