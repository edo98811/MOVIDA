# Module for Overview Dashboard with Grid Rendering
mod_dashboard_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Dashboard Module"),
    uiOutput(ns("grid_ui")) # Dynamic UI for the grid
  )
}

mod_dashboard_server <- function(id, elements_list) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output$grid_ui <- renderUI({
      if (is.null(elements_list) || length(elements_list) == 0) {
        return(h4("Add elements to the dashboard to see them here"))
      }
      
      # Render elements in a grid layout
      fluidRow(
        lapply(elements_list, function(element) {
          column(3, element) # Adjust column width as needed
        })
      )
    })
  })
}