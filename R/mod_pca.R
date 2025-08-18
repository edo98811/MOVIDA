#' module_plotting_lineplot UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_module_plotting_lineplot_ui <- function(id) {
  ns <- NS(id)
  tagList(
 
  )
}
    
#' module_plotting_lineplot Server Functions
#'
#' @noRd 
mod_module_plotting_lineplot_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_module_plotting_lineplot_ui("module_plotting_lineplot_1")
    
## To be copied in the server
# mod_module_plotting_lineplot_server("module_plotting_lineplot_1")
