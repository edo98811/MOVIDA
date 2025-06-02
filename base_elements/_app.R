library(shiny)

# Source UI and server scripts
source("_ui.R")
source("_server.R")

# Run the app

# Initialize once, probably in global.R or server.R at start
movida_data <- MovidaModelSingleton$get_instance(movida_list)

shinyApp(ui = ui, server = server)