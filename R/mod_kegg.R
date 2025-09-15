#' module_kegg UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_kegg_ui <- function(id) {
  ns <- NS(id)
  tagList(
 
  )
}
    
#' module_kegg Server Functions
#'
#' @noRd 
mod_kegg_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_module_kegg_ui("module_kegg_1")
    
## To be copied in the server
# mod_module_kegg_server("module_kegg_1")
