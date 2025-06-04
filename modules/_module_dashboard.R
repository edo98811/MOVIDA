source("modules/module_plot_container.R")


# Module for Overview Dashboard with Grid Rendering
mod_dashboard_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Dashboard Module"),
    uiOutput(ns("grid_ui")) # Dynamic UI for the grid
  )
}

mod_dashboard_server <- function(id, dashboard_elements, movida_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Render the dynamic UI for the dashboard grid
    output$grid_ui <- renderUI({
      elements <- dashboard_elements$elements

      if (length(elements) == 0) {
        return(h4("Add elements to the dashboard to see them here"))
      }

      do.call(layout_column_wrap, c(
        list(width = 1 / 3, heights_equal = "row"),
        lapply(names(elements), function(el_name) {
          list(mod_plot_ui(ns(el_name)))
        })
      ))
    })

    # Dynamically call module servers whenever dashboard_elements changes
    observeEvent(dashboard_elements$elements, {
      elements <- dashboard_elements$elements

      lapply(names(elements), function(el_name) {
        el_fun <- elements[[el_name]]
        message("created dashboard module server for: ", el_name)
        mod_plot_server(
          id = el_name,
          main_plotting_function = el_fun,
          dashboard_elements = dashboard_elements
        )
      })
    })
  })
}
