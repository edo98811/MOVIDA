library(shiny)
library(bslib)

server <- function(input, output, session, movida_data) {

  dashboard_elements <- reactiveValues(
    elements = list()  # Initialize an empty list to store dashboard elements
  )
  
  selected_row_source <- reactiveValues(
    selected = NULL,
    source = NULL,
    contrast = NULL
  )

  # Reactive data frame to store and update data dynamically
  reactive_data <- reactiveVal(data.frame())

  output$selected_feature <- renderUI({
    if (!is.null(selected_row_source$selected)) {
      div(
        style = "background-color: rgba(0, 128, 0, 0.8); height: 60px; display: flex; align-items: center; justify-content: space-between; padding: 20px 20px; color: white;",
        span(
          paste("Selected feature:", selected_row_source$selected, "from", selected_row_source$source)
        ),
        div(
          style = "margin-left: auto;",  # Push the button to the right
          actionButton("add_to_bookmarks", "Add to Bookmarks", style = "background-color: white; color: green; border: none;")
        ),
        hr()
      )
    } else {
      NULL  # Return NULL if selected is NULL
    }
  })

  # Call module servers
  mod_overview_server("overview", dashboard_elements, selected_row_source, movida_data)
  mod_sidebar_overview_server("sidebar_overview", dashboard_elements, selected_row_source, movida_data)
  mod_dashboard_server("dashboard", dashboard_elements, movida_data)  # Pass NULL or a list of elements as needed

  sidebar_content <- reactiveVal(FALSE)  # start with sidebar info hidden (set to FALSE)

  observeEvent(selected_row_source$contrast, {
    sidebar_content(FALSE)  # Hide sidebar content when contrast changes (resetting the selection)
  })

  observeEvent(selected_row_source$selected, {
    sidebar_content(TRUE)  # Hide sidebar content when contrast changes (resetting the selection)
  })

  # Render the sidebar UI based on the sidebar_content reactive value
  observeEvent(sidebar_content(), {
    output$sidebar_ui <- renderUI({
      if (sidebar_content()) {
        mod_sidebar_overview_ui("sidebar_overview")
      } else {
        # If sidebar is not visible, return an empty tagList
        tagList(h4("Select a row to show information"))
      }
    })
  })
}