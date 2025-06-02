library(shiny)

# Source UI and server scripts
source("_ui.R")
source("_server.R")

# Run the app
shinyApp(ui = ui, server = server)