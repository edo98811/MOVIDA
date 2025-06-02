mod_sidebar_overview_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Sidebar Overview - Main Panel Placeholder")
    # Add main panel UI elements here
    div(
      class = "base_container",
      plotOutput(ns("heatmap_plot_1")),
      plotOutput(ns("heatmap_plot_2"))
    )
  )
}

mod_sidebar_overview_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Initialize gargoyle for this session
    init("heatmap_event")
    
    # Reactive values to store event parameters
    event_params <- reactiveVal(NULL)
    
    # Watch the Gargoyle event and capture parameters
    watch("heatmap_event", {
      args <- list(...)
      event_params(args)
    }, session = session)
    
    # Create heatmap plots based on event parameters
    output$heatmap_plot_1 <- renderPlot({
      params <- event_params()
      if (!is.null(params)) {
        heatmap1 <- HeatmapBase$new(params$heatmap1)
        heatmap1$plot()
      }
    })
    
    output$heatmap_plot_2 <- renderPlot({
      params <- event_params()
      if (!is.null(params)) {
        heatmap2 <- HeatmapBase$new(params$heatmap2)
        heatmap2$plot()
      }
    })
  })
  })
}

# library(shiny)
# library(gargoyle)

# # UI function (empty or minimal, you can customize)
# myModuleUI <- function(id) {
#   ns <- NS(id)
#   tagList(
#     # Add UI elements here if needed
#     verbatimTextOutput(ns("event_info"))
#   )
# }

# # Server function
# myModuleServer <- function(id) {
#   moduleServer(id, function(input, output, session) {
#     # Initialize gargoyle for this session
#     init("my_event")
    
#     # Reactive values to store event parameters
#     event_params <- reactiveVal(NULL)
    
#     # Watch the Gargoyle event and capture parameters
#     watch("my_event", {
#       args <- list(...)
#       event_params(args)
#       # Your logic goes here
#       # For example:
#       # if (args$source == "section_A") { ... }
#     }, session = session)
    
#     # Example output to show event parameters (for testing)
#     output$event_info <- renderPrint({
#       event_params()
#     })
#   })
# }
