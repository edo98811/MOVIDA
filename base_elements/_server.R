library(shiny)
library(bslib)
library(gargoyle)

movida_data <- MovidaModelSingleton$get_instance(movida_list)

server <- function(input, output, session) {

  # Call module servers
  mod_overview_server("overview")
  mod_sidebar_overview_server("sidebar_overview")
  mod_dashboard_server("dashboard", elements_list = list())  # Pass NULL or a list of elements as needed
  
  # Init gargoyle events for this session
  init("toggle_sidebar")

  sidebar_content <- reactiveVal(TRUE)  # start visible (set to TRUE)

  # Listen for toggle_sidebar event to flip sidebar visibility
  on("toggle_sidebar", {
    sidebar_content(!sidebar_content())  # toggle the value
  })

  output$sidebar_ui <- renderUI({
    if (sidebar_content()) {
      mod_sidebar_overview_ui("sidebar_overview")
    } else {
      # If sidebar is not visible, return an empty tagList
      tagList(h4("Press a button to show information") )
    }
  })
}