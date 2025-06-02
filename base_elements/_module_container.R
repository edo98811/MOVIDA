library(shiny)

# UI function
mod_plot_ui <- function(id, dashboard = FALSE) {
  ns <- NS(id)
  layout_columns(
    column(
      width = 12,
      uiOutput(ns("dashboard_buttons"))
      if (!is.null(self$content$export_data())) {
        btns <- append(btns, actionButton(ns("btn_export_data"), "Export Data"))
      }
      if (!is.null(self$content$export_plot())) {
        btns <- append(btns, actionButton(ns("btn_export_plot"), "Export Plot"))
      }
    ),
    column(
      width = 12,
      plotOutput(ns("plot"))
    )
  )
}

# Server function
mod_plot_server <- function(id, r6_instance, r6_classes_list = NULL) {
  moduleServer(id, function(input, output, session) {
    output$plot <- renderPlot({
      req(r6_instance)
      # Call the render_plot method from the R6 instance
      r6_instance$render_plot()
    })

    output$dashboard_buttons <- renderUI({
      req(r6_instance)
      if (!is.null(r6_classes_list) && !is.null(r6_classes_list[[id]])) {
      actionButton(session$ns("btn_remove_dashboard"), "Remove from Dashboard")
      } else {
      actionButton(session$ns("btn_dashboard"), "Add to Dashboard")
      }
    })
    
    observeEvent(input$btn_dashboard, {
      req(r6_instance)
      if (!is.null(r6_classes_list)) {
        r6_classes_list[[id]] <- r6_instance
      }
    })
    
    observeEvent(input$btn_remove_dashboard, {
      req(r6_classes_list)
      r6_classes_list[[id]] <- NULL
    })
  })
}
